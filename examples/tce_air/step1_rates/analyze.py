"""Step 1 analysis: equilibrium rates of the O2 + N2 test case.

Produces, in results/step1_rates/:
  fig_progression.pdf     treatments of electronic energy; electronically
                            driven channel
  fig_eq_discretevib.pdf  discrete vibration
  table_reservoir.tex/.txt  equilibrium rows and frozen-electronic rows
  numbers.json              every number the paper quotes from this step
  results.txt               PASS/FAIL checks against independent references
"""
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
from runs import runs, SIZES, STEP  # noqa: E402
from pkg import size_params  # noqa: E402

V = 1.0e-12                 # box volume (m^3)
NRHO = 7.07043e22           # total number density (m^-3)
DT = 1.0e-9


def rate(run):
    """rate coefficient from the tally of one frozen reservoir run"""
    n = pkg.tallies(run.log).get("O2 + N2 --> O + O + N2", 0)
    NP, NSTEP = run.variables["NP"], run.variables["NSTEP"]
    fnum = NRHO * V / NP
    n1 = n2 = 0.5 * NRHO
    return n * fnum / (n1 * n2 * V * NSTEP * DT), n


def collect(rl, group, treat, T, TE=None, rot=False):
    sel = [r for r in rl if r.tag.startswith("%s_%s_T%d%s%s_s"
                                             % (group, treat, T,
                                                "_TE%d" % TE if TE else "",
                                                "_rot" if rot else ""))]
    sel = [r for r in sel if r.done()]
    if not sel:
        return None
    ks, ns = zip(*[rate(r) for r in sel])
    m, sd = pkg.mean_sd(list(ks))
    sem = sd / math.sqrt(len(ks)) if len(ks) > 1 else 0.0
    poisson = m / math.sqrt(max(sum(ns), 1)) if sum(ns) > 0 else 0.0
    return dict(k=m, err=max(sem, poisson), n=sum(ns), nseed=len(ks))


def analyze(size="quick"):
    rep = pkg.Report(STEP + " (" + size + ")")
    out = os.path.join(pkg.results_dir(), STEP)
    os.makedirs(out, exist_ok=True)
    rl = runs(size)
    D = R.Data(elec=os.path.join(HERE, "o2n2.elec"),
               tce=os.path.join(HERE, "o2n2.tce"))
    rx = D.rx[0]
    O2, N2 = D.elec["O2"], D.elec["N2"]
    z0 = D.zrot("O2", "N2")
    om = D.pair("O2", "N2")[1]
    eta, Ea, A = rx["eta"], rx["Ea"], rx["A"]
    ex = lambda lad: [0.0 if t == 0 else 1.0 for t, g in lad]
    # refit of the per-state DOF model, the one the deck reads
    Rf = R.Data(elec=os.path.join(HERE, "o2n2.elec"),
                tce=os.path.join(HERE, "o2n2_refit.tce")).rx[0]
    num = {}

    def kref(T):
        return R.arrhenius(rx, T)

    def Rrefit(T):
        return (R.R_inclusion(T, O2, N2, z0, om, Rf["eta"], Rf["Ea"],
                              ex(O2), ex(N2)) * R.arrhenius(Rf, T) / kref(T))

    curves = {
        "incl": lambda T: R.R_inclusion(T, O2, N2, z0, om, eta, Ea),
        "dof1": lambda T: R.R_inclusion(T, O2, N2, z0, om, eta, Ea,
                                        ex(O2), ex(N2)),
        "dof1refit": Rrefit,
        "micro": lambda T: 1.0,
    }

    # ---------------- Fig. 1 left: electronic treatments, vibration left out
    rep.info("Fig. 1: k/k_Arr with vibration left out (z = %g, omega = %.3f)"
             % (z0, om))
    fig1 = {}
    for tr in ("micro", "incl", "dof1", "dof1refit"):
        for T in (5000, 10000, 20000, 40000):
            c = collect(rl, "fig1", tr, T)
            if c is None:
                continue
            ref = curves[tr](T)
            r, e = c["k"] / kref(T), c["err"] / kref(T)
            fig1[(tr, T)] = (r, e)
            num["fig1_%s_%d" % (tr, T)] = [r, 2 * e, ref]
            rep.check(abs(r - ref) <= 3 * e + 0.01 * ref,
                      "%-9s T=%5d K: measured %.3f +/- %.3f (2 sd), "
                      "quadrature %.3f" % (tr, T, r, 2 * e, ref),
                      meas=r, ref=ref, sd=e, tol=0.01 * ref,
                      label="no vibration, %s, %d K" % (tr, T))

    # ---------------- Fig. 1 right: electronically driven channel
    rep.info("Fig. 1 right: Ttr = 2000 K, electronic populations frozen at "
             "40 000 K (reaction counts per seed)")
    T, TE = 2000.0, 40000.0
    pred = {"micro": R.micro_ratio(T, z0, om, eta, Ea, [], [(O2, TE), (N2, TE)]),
            "incl": R.R_inclusion(T, O2, N2, z0, om, eta, Ea, Tpop=TE),
            "dof1": R.R_inclusion(T, O2, N2, z0, om, eta, Ea, ex(O2), ex(N2),
                                  Tpop=TE),
            "excl": 1.0}
    chan = {}
    NPs, NST = (size_params(STEP, SIZES, size)[k] for k in ("NP", "NSTEP"))
    fnum = NRHO * V / NPs
    count_per_k = (0.5 * NRHO) ** 2 * V * NST * DT / fnum
    for tr in ("excl", "dof1", "micro", "incl"):
        c = collect(rl, "chan", tr, 2000, TE=40000)
        if c is None:
            continue
        exp_n = pred[tr] * kref(T) * count_per_k
        meas_n = c["n"] / c["nseed"]
        chan[tr] = (meas_n, exp_n)
        num["chan_%s" % tr] = [meas_n, exp_n]
        if tr == "excl":
            rep.check(c["n"] < 5, "electronic energy left out: %d reactions "
                      "(expected %.2g per seed)" % (c["n"], exp_n))
        else:
            rep.check(abs(meas_n - exp_n) <= 3 * math.sqrt(exp_n / c["nseed"])
                      + 0.01 * exp_n,
                      "%-5s %.0f reactions per seed, predicted %.0f (%+.2f%%)"
                      % (tr, meas_n, exp_n, 100 * (meas_n / exp_n - 1)),
                      meas=meas_n, ref=exp_n,
                      sd=math.sqrt(exp_n / c["nseed"]), tol=0.01 * exp_n,
                      label="2000 K, T_el 40 000 K channel, %s" % tr)

    # ---------------- Table VI top and Fig. 3: discrete vibration
    rep.info("Table VI: k/k_Arr with discrete vibration")
    th = (D.sp["O2"]["vibtemp"], D.sp["N2"]["vibtemp"])
    tab6 = {}
    for T in (10000, 15000, 20000):
        inst_q = R.R_instdof(T, th[0], th[1], z0, om, eta, Ea)
        num["instdof_quad_%d" % T] = inst_q
        for tr in ("inst", "incl", "excl", "microel", "micro", "vibonly"):
            c = collect(rl, "tab6", tr, T)
            if c is None:
                continue
            tab6[(tr, T)] = (c["k"] / kref(T), c["err"] / kref(T))
            num["tab6_%s_%d" % (tr, T)] = [c["k"] / kref(T),
                                           2 * c["err"] / kref(T)]
        for tr, ref, what in (("inst", inst_q, "instantaneous DOF vs quadrature"),
                              ("excl", inst_q, "electronic energy left out vs "
                                               "instantaneous-DOF quadrature"),
                              ("micro", 1.0, "joint vib+elec factor vs input"),
                              ("vibonly", 1.0, "vibrational ladder factor vs "
                                               "input")):
            if (tr, T) in tab6:
                r, e = tab6[(tr, T)]
                rep.check(abs(r - ref) <= 3 * e + 0.01 * ref,
                          "T=%5d K %s: %.3f +/- %.3f, expected %.3f"
                          % (T, what, r, 2 * e, ref),
                          meas=r, ref=ref, sd=e, tol=0.01 * ref)
        if ("microel", T) in tab6 and ("inst", T) in tab6:
            r = tab6[("microel", T)][0] / tab6[("inst", T)][0]
            rep.info("T=%5d K microcanonical electronic factor alone / "
                     "instantaneous-DOF baseline = %.3f" % (T, r))
    c = collect(rl, "tab6", "micro", 15000, rot=True)
    if c:
        r, e = c["k"] / kref(15000), c["err"] / kref(15000)
        num["tab6_micro_rot_15000"] = [r, 2 * e]
        rep.check(abs(r - 1) <= 3 * e + 0.01, "discrete rotation, joint "
                  "factor at 15 kK: %.3f +/- %.3f, expected 1" % (r, 2 * e),
                  meas=r, ref=1.0, sd=e, tol=0.01)

    # ---------------- Table VI bottom: frozen electronic populations
    rep.info("Table VI bottom: Ttr = 15 kK, electronic populations frozen, "
             "relative to the baseline equilibrium rate (inst, 15 kK)")
    base = collect(rl, "tab6", "inst", 15000)
    frz = {}
    for TE in (100, 40000):
        frz[("analytic", TE)] = R.micro_ratio(15000.0, z0, om, eta, Ea, [],
                                              [(O2, float(TE)), (N2, float(TE))])
        for tr in ("microel", "incl", "excl"):
            c = collect(rl, "frz", tr, 15000, TE=TE)
            if c is None or base is None:
                continue
            r = c["k"] / base["k"]
            e = r * math.hypot(c["err"] / c["k"], base["err"] / base["k"])
            frz[(tr, TE)] = (r, e)
            num["frz_%s_%d" % (tr, TE)] = [r, 2 * e]
        num["frz_analytic_%d" % TE] = frz[("analytic", TE)]
    for TE in (100, 40000):
        if ("excl", TE) in frz:
            r, e = frz[("excl", TE)]
            rep.check(abs(r - 1) <= 3 * e + 0.01, "TE=%5d K electronic energy "
                      "left out: %.3f +/- %.3f, expected 1" % (TE, r, 2 * e),
                      meas=r, ref=1.0, sd=e, tol=0.01)
    if ("incl", 100) in frz:
        r, e = frz[("incl", 100)]
        rep.check(abs(r - 1) <= 3 * e + 0.01, "TE=100 K inclusion (ground "
                  "states only): %.3f +/- %.3f, expected 1" % (r, 2 * e),
                  meas=r, ref=1.0, sd=e, tol=0.01)

    _write_table6(out, tab6, frz)
    _fig1(out, fig1, chan, curves)
    _fig3(out, tab6, th, z0, om, eta, Ea)
    with open(os.path.join(out, "numbers.json"), "w") as fh:
        json.dump({k: v for k, v in num.items()}, fh, indent=1)
    return rep.write()


def _fmt(d, key):
    return "%.3f" % d[key][0] if key in d else "--"


def _write_table6(out, tab6, frz):
    rows = [("instantaneous DOF, no electronic coupling", "inst"),
            ("\\quad electronic energy included", "incl"),
            ("\\quad electronic energy excluded", "excl"),
            ("\\quad microcanonical factor, electronic ladder only", "microel"),
            ("joint vibrational--electronic factor", "micro")]
    L = ["equilibrium, $k/k_{\\mathrm{Arr}}$ & 10~kK & 15~kK & 20~kK \\\\",
         "\\hline"]
    for name, tr in rows:
        L.append("%s & %s \\\\" % (name, " & ".join(_fmt(tab6, (tr, T))
                                                    for T in (10000, 15000,
                                                              20000))))
    L += ["\\hline",
          "frozen $\\Tel$, $k/k_{\\mathrm{baseline}}$ & $\\Tel = 100$~K & & "
          "$\\Tel = 40\\,000$~K \\\\", "\\hline",
          "analytic, constant $z$ & %.3f & & %.3f \\\\"
          % (frz[("analytic", 100)], frz[("analytic", 40000)])]
    for name, tr in (("microcanonical factor, electronic ladder only",
                      "microel"),
                     ("electronic energy included", "incl"),
                     ("electronic energy excluded", "excl")):
        L.append("%s & %s & & %s \\\\" % (name, _fmt(frz, (tr, 100)),
                                          _fmt(frz, (tr, 40000))))
    open(os.path.join(out, "table_reservoir.tex"), "w").write("\n".join(L) + "\n")
    open(os.path.join(out, "table_reservoir.txt"), "w").write(
        "\n".join(l.replace("\\\\", "").replace("&", "|") for l in L) + "\n")


def _fig1(out, fig1, chan, curves):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import paper_style as ps
    ps.apply()
    fig, (a, b) = plt.subplots(1, 2, figsize=ps.DOUBLE,
                               gridspec_kw=dict(width_ratios=[1.6, 1]))
    Tc = np.linspace(4000.0, 45000.0, 120)
    sty = {"incl": (ps.VERM, "o", "electronic energy included"),
           "dof1": (ps.ORANGE, "s", "one DOF per excited state"),
           "dof1refit": (ps.PURPLE, "D", "same, Arrhenius refit"),
           "micro": (ps.BLUE, "^", "microcanonical factor")}
    for tr, (col, mk, lab) in sty.items():
        a.plot(Tc / 1e3, [curves[tr](T) for T in Tc], color=col, lw=1.4,
               label=lab)
        pts = [(T, *fig1[(tr, T)]) for T in (5000, 10000, 20000, 40000)
               if (tr, T) in fig1]
        if pts:
            T_, r_, e_ = zip(*pts)
            a.errorbar(np.array(T_) / 1e3, r_, yerr=2 * np.array(e_), fmt=mk,
                       color=col, mfc="white", ms=4.5, capsize=2, lw=0)
    a.axhline(1.0, color=ps.GRAY, lw=0.8, ls=":")
    a.set_xlabel("$T$ (kK)")
    a.set_ylabel("$k_{\\mathrm{model}}/k$")
    a.set_xlim(4, 45)
    a.legend(loc="upper right", fontsize=7)
    ps.panel(a, "(a)")
    names = [("excl", "left out"), ("dof1", "1 DOF\nper state"),
             ("micro", "micro-\ncanonical"), ("incl", "included")]
    names = [(k, n) for k, n in names if k in chan]
    x = np.arange(len(names))
    b.bar(x - 0.18, [chan[k][0] / 1e3 for k, _ in names], 0.36,
          color=ps.BLUE, label="SPARTA")
    b.bar(x + 0.18, [chan[k][1] / 1e3 for k, _ in names], 0.36,
          color=ps.GRAY, alpha=0.5, label="quadrature")
    b.set_xticks(x)
    b.set_xticklabels([n for _, n in names], fontsize=7.5)
    b.set_ylabel("reactions per run ($10^3$)")
    b.legend(loc="upper left", fontsize=7)
    b.grid(axis="x", visible=False)
    ps.panel(b, "(b)", "ne")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig_progression.pdf"))
    plt.close(fig)


def _fig3(out, tab6, th, z0, om, eta, Ea):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import paper_style as ps
    ps.apply()
    fig, a = plt.subplots(figsize=ps.SINGLE)
    Tc = np.linspace(8000.0, 25000.0, 35)
    a.plot(Tc / 1e3, [R.R_instdof(T, th[0], th[1], z0, om, eta, Ea)
                      for T in Tc], color=ps.VERM, lw=1.4,
           label="instantaneous DOF, quadrature")
    for tr, col, mk, lab in (("inst", ps.VERM, "s", "instantaneous DOF"),
                             ("microel", ps.ORANGE, "D",
                              "microcanonical, electronic ladder only"),
                             ("micro", ps.BLUE, "o",
                              "joint vibrational-electronic factor")):
        pts = [(T, *tab6[(tr, T)]) for T in (10000, 15000, 20000)
               if (tr, T) in tab6]
        if pts:
            T_, r_, e_ = zip(*pts)
            a.errorbar(np.array(T_) / 1e3, r_, yerr=2 * np.array(e_), fmt=mk,
                       color=col, mfc="white", ms=5, capsize=2, lw=0,
                       label=lab)
    a.axhline(1.0, color=ps.GRAY, lw=0.8, ls=":")
    a.set_xlabel("$T$ (kK)")
    a.set_ylabel("$k/k_{\\mathrm{Arr}}$")
    a.legend(loc="upper right", fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig_eq_discretevib.pdf"))
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(analyze(sys.argv[1] if len(sys.argv) > 1 else "quick"))
