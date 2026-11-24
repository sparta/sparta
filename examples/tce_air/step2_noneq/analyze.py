"""Step 2 analysis: non-equilibrium response (Sec. VII B).

Produces, in results/step2_noneq/:
  fig_noneq.pdf        rate relative to the input rate, (a)
                         electronic populations frozen at TE, (b) vibrational
                         populations frozen at TV; one normalization
  fig_reservoir.pdf    (a) equilibrium checks of the factor (from
                         step 1), (b) measured against analytic frozen-ladder
                         factors
  fig_qk_compare.pdf   microcanonical factor, quantum-kinetic model
                         and Park's two-temperature forms
  table_qct.tex/.txt    model against QCT (Bender et al. 2015)
  numbers.json, results.txt
"""
import csv
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)
import pkg                 # noqa: E402
import tce_ref as R        # noqa: E402
sys.modules.pop("runs", None)  # this step's runs.py, not another step's
from runs import runs, STEP, T, TVS, TES  # noqa: E402

V, NRHO, DT = 1.0e-12, 7.07043e22, 1.0e-9


def rate(run):
    n = pkg.tallies(run.log).get("O2 + N2 --> O + O + N2", 0)
    NP, NSTEP = run.variables["NP"], run.variables["NSTEP"]
    fnum = NRHO * V / NP
    return n * fnum / ((0.5 * NRHO) ** 2 * V * NSTEP * DT), n


def collect(rl, prefix):
    sel = [r for r in rl if r.tag.startswith(prefix + "_s") and r.done()]
    if not sel:
        return None
    ks, ns = zip(*[rate(r) for r in sel])
    m, sd = pkg.mean_sd(list(ks))
    sem = sd / math.sqrt(len(ks)) if len(ks) > 1 else 0.0
    return dict(k=m, err=max(sem, m / math.sqrt(max(sum(ns), 1))),
                n=sum(ns))


def analyze(size="quick"):
    rep = pkg.Report(STEP + " (" + size + ")")
    out = os.path.join(pkg.results_dir(), STEP)
    os.makedirs(out, exist_ok=True)
    rl = runs(size)
    D = R.Data(elec=os.path.join(HERE, "o2n2.elec"),
               tce=os.path.join(HERE, "o2n2.tce"))
    rx = D.rx[0]
    z0, om = D.zrot("O2", "N2"), D.pair("O2", "N2")[1]
    eta, Ea_K = rx["eta"], rx["Ea"] / R.KB
    th = (D.sp["O2"]["vibtemp"], D.sp["N2"]["vibtemp"])
    J = R.JointLadder(z0, om, eta, Ea_K, thetas=th,
                      elecs=(D.elec["O2"], D.elec["N2"]))
    Jv = R.JointLadder(z0, om, eta, Ea_K, thetas=th, du=20.0)
    ka = R.arrhenius(rx, T)
    num = {}

    # ---------------- frozen vibration
    vib = {}
    for TV in TVS:
        ref = J.rate(T, Tv=TV)
        c = collect(rl, "vib_micro_TV%d" % TV)
        if c:
            r, e = c["k"] / ka, c["err"] / ka
            vib[TV] = (r, e, ref)
            num["vib_micro_%d" % TV] = [r, 2 * e, ref]
            rep.check(abs(r - ref) <= 3 * e + 0.01 * ref,
                      "TV=%5d K microcanonical: %.3f +/- %.3f, analytic %.3f"
                      % (TV, r, 2 * e, ref),
                      meas=r, ref=ref, sd=e, tol=0.01 * ref)
        refi = R.R_instdof(T, th[0], th[1], z0, om, eta, rx["Ea"], Tvib=TV)
        c = collect(rl, "vib_inst_TV%d" % TV)
        if c:
            r, e = c["k"] / ka, c["err"] / ka
            vib[("inst", TV)] = (r, e, refi)
            num["vib_inst_%d" % TV] = [r, 2 * e, refi]
            rep.check(abs(r - refi) <= 3 * e + 0.01 * refi,
                      "TV=%5d K instantaneous DOF: %.3f +/- %.3f, quadrature "
                      "%.3f" % (TV, r, 2 * e, refi),
                      meas=r, ref=refi, sd=e, tol=0.01 * refi)

    # ---------------- frozen electronic
    el = {}
    for TE in TES:
        ref = J.rate(T, Te=TE)
        c = collect(rl, "el_micro_TE%d" % TE)
        if c:
            r, e = c["k"] / ka, c["err"] / ka
            el[TE] = (r, e, ref)
            num["el_micro_%d" % TE] = [r, 2 * e, ref]
            rep.check(abs(r - ref) <= 3 * e + 0.01 * ref,
                      "TE=%5d K microcanonical: %.3f +/- %.3f, analytic %.3f"
                      % (TE, r, 2 * e, ref),
                      meas=r, ref=ref, sd=e, tol=0.01 * ref)
        c = collect(rl, "el_incl_TE%d" % TE)
        if c:
            el[("incl", TE)] = (c["k"] / ka, c["err"] / ka)
            num["el_incl_%d" % TE] = [c["k"] / ka, 2 * c["err"] / ka]

    # ---------------- QK and Park (calculations)
    qk = lambda Tt, Tv: R.qk_over_kcoll(Tt, Tv, th[0], Ea_K, om)
    for Tt in (10000.0, 15000.0, 20000.0, 30000.0):
        num["qk_over_park_%d" % Tt] = (qk(Tt, Tt) * D.kcoll("O2", "N2", Tt)
                                      / R.arrhenius(rx, Tt))
    num["qk_mult_Tv3000"] = qk(T, 3000.0) / qk(T, T)
    num["micro_mult_Tv3000"] = J.rate(T, Tv=3000.0)
    num["micro_vibonly_mult_Tv3000"] = Jv.rate(T, Tv=3000.0)
    num["park_q05_Tv3000"] = R.park_ratio(T, 3000.0, 0.5, eta, Ea_K)
    rep.info("QK/Park fit at 10/15/20/30 kK: %s" % ", ".join(
        "%.2f" % num["qk_over_park_%d" % t] for t in (10000, 15000, 20000,
                                                      30000)))
    rep.info("at TV = 3000 K: QK x%.2f, microcanonical x%.3f (vibrational "
             "ladder alone x%.3f), Park q=0.5 x%.4f"
             % (num["qk_mult_Tv3000"], num["micro_mult_Tv3000"],
                num["micro_vibonly_mult_Tv3000"], num["park_q05_Tv3000"]))

    # ---------------- QCT table (calculation)
    qct = _qct_table(rep, num)

    _fig4(out, vib, el, J, th, z0, om, eta, rx["Ea"])
    _fig5(out, vib, el, size)
    _fig6(out, vib, J, Jv, qk, eta, Ea_K)
    _table5(out, qct)
    with open(os.path.join(out, "numbers.json"), "w") as fh:
        json.dump(num, fh, indent=1, default=float)
    return rep.write()


def _qct_table(rep, num):
    Dair = R.Data()
    rx = Dair.reaction("N2 + N2 --> N + N + N2")
    z0, om = Dair.zrot("N2", "N2"), Dair.pair("N2", "N2")[1]
    th = Dair.sp["N2"]["vibtemp"]
    J = R.JointLadder(z0, om, rx["eta"], rx["Ea"] / R.KB, thetas=(th, th),
                      du=th / 16.0)
    qct = {}
    for row in csv.DictReader(l for l in open(os.path.join(
            HERE, "bender2015_k_T_Tv.csv")) if not l.startswith("#")):
        Tt, Tv = float(row["T_K"]), float(row["Tv_K"])
        qct[(Tt, Tv)] = float(row["ratio_to_thermal"])
    tab = {}
    for (Tt, Tv), q in qct.items():
        m = J.rate(Tt, Tv=Tv) / J.rate(Tt)
        tab[(Tt, Tv)] = (m, q)
    shown = (10000.0, 13000.0, 20000.0)          # the columns of Table V
    frac_cold = [math.log(m) / math.log(q) for (Tt, Tv), (m, q) in tab.items()
                 if Tt in shown and Tv < Tt and q > 0]
    frac_hot = [math.log(m) / math.log(q) for (Tt, Tv), (m, q) in tab.items()
                if Tt in shown and Tv > Tt and q > 0]
    num["qct_logfrac_cold"] = [min(frac_cold), max(frac_cold)]
    num["qct_logfrac_hot"] = [min(frac_hot), max(frac_hot)]
    num["qct_table"] = {"%d_%d" % k: v for k, v in tab.items()}
    rep.info("model / QCT logarithmic response: vibrationally cold %.2f-%.2f, "
             "hot %.2f-%.2f (N2 + N2 fit of data/air.tce: eta %g, Ea %.4e J)"
             % (min(frac_cold), max(frac_cold), min(frac_hot), max(frac_hot),
                rx["eta"], rx["Ea"]))
    return tab


def _table5(out, tab):
    Ts = sorted({k[0] for k in tab})
    Tvs = sorted({k[1] for k in tab})
    show_T = [t for t in (10000.0, 13000.0, 20000.0) if t in Ts]
    L = ["$\\Tvib$ (K) & " + " & ".join("model & QCT" for _ in show_T)
         + " \\\\", "\\hline"]
    for Tv in Tvs:
        cells = []
        for Tt in show_T:
            if (Tt, Tv) in tab:
                m, q = tab[(Tt, Tv)]
                cells += ["%.2g" % m, "%.2g" % q]
            else:
                cells += ["--", "--"]
        L.append("%g & %s \\\\" % (Tv, " & ".join(cells)))
    open(os.path.join(out, "table_qct.tex"), "w").write("\n".join(L) + "\n")
    open(os.path.join(out, "table_qct.txt"), "w").write(
        "\n".join(l.replace("\\\\", "").replace("&", "|") for l in L) + "\n")


def _setup():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import paper_style as ps
    ps.apply()
    return plt, ps


def _fig4(out, vib, el, J, th, z0, om, eta, Ea):
    plt, ps = _setup()
    fig, (a, b) = plt.subplots(1, 2, figsize=ps.DOUBLE, sharey=False)
    Tc = np.geomspace(100.0, 45000.0, 24)
    a.plot(Tc / 1e3, [J.rate(T, Te=t) for t in Tc], color=ps.BLUE,
           label="microcanonical, analytic")
    _pts(a, [(t, *el[t][:2]) for t in TES if t in el], ps.BLUE, "o",
         "microcanonical, SPARTA")
    _pts(a, [(t, *el[("incl", t)]) for t in TES if ("incl", t) in el],
         ps.VERM, "s", "inclusion + instantaneous DOF, SPARTA")
    a.set_xscale("log")
    a.set_xlabel("$T_{\\mathrm{elec}}$ (kK)")
    a.set_ylabel("$k/k_{\\mathrm{Arr}}(T)$")
    a.legend(fontsize=6.8, loc="upper left")
    ps.panel(a, "(a)", "se")
    Tc = np.geomspace(300.0, 45000.0, 24)
    b.plot(Tc / 1e3, [J.rate(T, Tv=t) for t in Tc], color=ps.BLUE,
           label="microcanonical, analytic")
    b.plot(Tc / 1e3, [R.R_instdof(T, th[0], th[1], z0, om, eta, Ea, Tvib=t)
                      for t in Tc], color=ps.VERM, ls="--",
           label="instantaneous DOF, quadrature")
    _pts(b, [(t, *vib[t][:2]) for t in TVS if t in vib], ps.BLUE, "o",
         "microcanonical, SPARTA")
    _pts(b, [(t, *vib[("inst", t)][:2]) for t in TVS if ("inst", t) in vib],
         ps.VERM, "s", "instantaneous DOF, SPARTA")
    b.set_xscale("log")
    b.set_xlabel("$T_{\\mathrm{vib}}$ (kK)")
    b.set_ylabel("$k/k_{\\mathrm{Arr}}(T)$")
    b.legend(fontsize=6.8, loc="upper left")
    ps.panel(b, "(b)", "se")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig_noneq.pdf"))
    plt.close(fig)


def _pts(ax, pts, col, mk, lab):
    if not pts:
        return
    x, y, e = zip(*pts)
    ax.errorbar(np.array(x) / 1e3, y, yerr=2 * np.array(e), fmt=mk, color=col,
                mfc="white", ms=4.5, capsize=2, lw=0, label=lab)


def _fig5(out, vib, el, size):
    """(a) equilibrium checks from step 1; (b) measured vs analytic"""
    plt, ps = _setup()
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "step1_runs", os.path.join(os.path.dirname(HERE), "step1_rates",
                                   "runs.py"))
    s1runs = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(s1runs)
    rl1 = s1runs.runs(size)

    class s1:      # step 1 reservoirs: same rate formula as this step
        @staticmethod
        def collect(rl, group, treat, Tt, rot=False):
            return collect(rl, "%s_%s_T%d%s" % (group, treat, Tt,
                                                "_rot" if rot else ""))
    D = R.Data(tce=os.path.join(HERE, "o2n2.tce"),
               elec=os.path.join(HERE, "o2n2.elec"))
    fig, (a, b) = plt.subplots(1, 2, figsize=ps.DOUBLE)
    for k, (tr, rot, lab, mk, col) in enumerate((
            ("vibonly", False, "vib", "o", ps.BLUE),
            ("micro", False, "vib+elec", "s", ps.GREEN),
            ("micro", True, "rot discrete", "D", ps.PURPLE))):
        pts = []
        for Tt in (10000, 15000, 20000):
            c = s1.collect(rl1, "tab6", tr, Tt, rot=rot)
            if c:
                ka = R.arrhenius(D.rx[0], Tt)
                pts.append((Tt + (k - 1) * 250, c["k"] / ka, c["err"] / ka))
        _pts(a, pts, col, mk, lab)
    a.axhline(1.0, color=ps.GRAY, lw=0.8, ls=":")
    a.set_xlabel("$T$ (kK)")
    a.set_ylabel("$k/k_{\\mathrm{Arr}}$")
    a.legend(fontsize=7)
    ps.panel(a, "(a)")
    pts_v = [(vib[t][2], vib[t][0], vib[t][1]) for t in TVS if t in vib]
    pts_e = [(el[t][2], el[t][0], el[t][1]) for t in TES if t in el]
    for pts, col, mk, lab in ((pts_v, ps.BLUE, "o", "vibration frozen"),
                              (pts_e, ps.ORANGE, "s", "electronic frozen")):
        if pts:
            x, y, e = zip(*pts)
            b.errorbar(x, y, yerr=2 * np.array(e), fmt=mk, color=col,
                       mfc="white", ms=4.5, capsize=2, lw=0, label=lab)
    lim = [0.4, 2.2]
    b.plot(lim, lim, color=ps.GRAY, lw=0.8, ls=":")
    b.set_xlim(lim)
    b.set_ylim(lim)
    b.set_xlabel("analytic factor")
    b.set_ylabel("SPARTA factor")
    b.legend(fontsize=7, loc="upper left")
    ps.panel(b, "(b)", "se")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig_reservoir.pdf"))
    plt.close(fig)


def _fig6(out, vib, J, Jv, qk, eta, Ea_K):
    plt, ps = _setup()
    fig, a = plt.subplots(figsize=ps.SINGLE)
    Tc = np.geomspace(300.0, 45000.0, 24)
    a.plot(Tc / 1e3, [J.rate(T, Tv=t) for t in Tc], color=ps.BLUE,
           label="microcanonical (vib+elec)")
    a.plot(Tc / 1e3, [Jv.rate(T, Tv=t) / Jv.rate(T) for t in Tc],
           color=ps.SKY, ls="-.", label="microcanonical, vibration only")
    a.plot(Tc / 1e3, [qk(T, t) / qk(T, T) for t in Tc], color=ps.GREEN,
           ls="--", label="quantum-kinetic (Bird)")
    for q, ls in ((0.5, ":"), (0.7, (0, (5, 2, 1, 2)))):
        a.plot(Tc / 1e3, [R.park_ratio(T, t, q, eta, Ea_K) for t in Tc],
               color=ps.VERM, ls=ls, label="Park, $q = %.1f$" % q)
    _pts(a, [(t, vib[t][0], vib[t][1]) for t in TVS if t in vib], ps.BLUE,
         "o", "SPARTA")
    a.set_xscale("log")
    a.set_yscale("log")
    a.set_ylim(1e-3, 5)
    a.set_xlabel("$T_{\\mathrm{vib}}$ (kK)")
    a.set_ylabel("$k(T, T_{\\mathrm{vib}})/k(T)$")
    a.legend(fontsize=6.8, loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig_qk_compare.pdf"))
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(analyze(sys.argv[1] if len(sys.argv) > 1 else "quick"))
