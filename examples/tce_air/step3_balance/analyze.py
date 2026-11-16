"""Step 3 analysis: detailed balance of the derived reverse reactions.

Produces, in results/step3_balance/:
  fig02_zeldovich.pdf   Fig. 7: (a) reverse rate of NO + N -> N2 + O
                        (derived, legacy fit, SPARTA); (b) K_eq implied by
                        the legacy forward/reverse fits relative to the
                        partition-function K_eq; (c) forward/reverse tally
                        ratio relative to K_eq
  table7_recomb.tex/.txt  Table VII (recombination configurations)
  fig_mechanism.pdf     the 28 dissociation and exchange pairs (and the 3
                        associative-ionization pairs) of the full mechanism
                        at 15 and 8 kK
  mechanism.txt, numbers.json, results.txt
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
from runs import runs, STEP, ZELD_T, REC, NLOW, NHIGH  # noqa: E402

V, DT = 1.0e-12, 1.0e-9


def sums(rl, prefix):
    sel = [r for r in rl if r.tag.startswith(prefix + "_s") and r.done()]
    if not sel:
        return None, []
    tot = {}
    for r in sel:
        for k, v in pkg.tallies(r.log).items():
            tot[k] = tot.get(k, 0) + v
    return tot, sel


def mean_T(sel):
    """translational temperature of the runs, averaged over time (stats
    rows after the first) and over the runs: the reference values are
    evaluated at it, since a closed reservoir settles a few tens of kelvin
    away from its nominal temperature"""
    vals = []
    for r in sel:
        head, rows = pkg.stats(r.log)
        i = head.index("c_temp")
        vals.append(np.mean([row[i] for row in rows[1:]]))
    return float(np.mean(vals))


RES_SPC = ("N2", "O", "NO", "N")
MECH_SPC = ("N2", "O2", "NO", "N", "O", "N2+", "O2+", "NO+", "N+", "O+", "e")


def composition(run, names):
    """mole fractions of the sample: the species counts of the first stats
    row (create_particles draws each particle's species at random, so they
    differ from the nominal fractions by about 1/sqrt(N); the reservoirs
    are frozen, so the composition does not change during the run)"""
    head, rows = pkg.stats(run.log)
    i = head.index("c_spc[1]")
    c = rows[0][i:i + len(names)]
    return dict(zip(names, [v / sum(c) for v in c]))


def run_T(run):
    """time-averaged translational temperature of one run"""
    head, rows = pkg.stats(run.log)
    i = head.index("c_temp")
    return float(np.mean([row[i] for row in rows[1:]]))


def ratio_err(F, B, per_seed):
    """ratio F/B with the larger of the Poisson and the seed error"""
    r = F / B
    e = r * math.sqrt(1.0 / F + 1.0 / B)
    if len(per_seed) > 1:
        m, sd = pkg.mean_sd(per_seed)
        e = max(e, r * sd / m / math.sqrt(len(per_seed)))
    return r, e


def check_lines(rep):
    """the reaction files of this step must agree with data/air.tce"""
    air = [" ".join(l.split()) for l in open(os.path.join(R.DATA, "air.tce"))]
    for f in ("zeldovich.tce", "recomb_atom.tce", "recomb_atom_fwd.tce",
              "recomb_mol.tce"):
        lines = [" ".join(l.split()) for l in open(os.path.join(HERE, f))
                 if l.strip() and not l.startswith("#")]
        ok = all(l in air for l in lines)
        rep.check(ok, "%s agrees with data/air.tce" % f)


def analyze(size="quick"):
    rep = pkg.Report(STEP + " (" + size + ")")
    out = os.path.join(pkg.results_dir(), STEP)
    os.makedirs(out, exist_ok=True)
    rl = runs(size)
    D = R.Data()
    num = {}
    check_lines(rep)

    # ---------------- Zeldovich pair
    fwd = D.reaction("N2 + O --> NO + N")
    # the derived reverse NO + N -> N2 + O is barrierless and its probability
    # exceeds one at low collision energy, where it is capped; the cap is
    # not carried into the forward table (Sec. IV D), so the reverse runs
    # slightly below k_f/K_eq and (f/b)/K_eq = 1/(capped fraction)
    XR = R.ExchangeReverse(D, fwd)
    zel = {}
    for T in ZELD_T:
        tot, sel = sums(rl, "zeld_T%d" % T)
        if not tot:
            continue
        Tn, T = T, mean_T(sel)
        F, B = tot.get(fwd["label"], 0), tot.get("NO + N --> N2 + O", 0)
        keq = D.keq(fwd, T)
        NP, NST = sel[0].variables["NP"], sel[0].variables["NSTEP"]
        fnum = NLOW * V / NP
        # each seed's sample: f/b = K_eq x_N2 x_O/(x_NO x_N), and the
        # reverse rate is its tally over n_NO n_N V t
        Bc, kden, per = 0.0, 0.0, []
        for run in sel:
            tl = pkg.tallies(run.log)
            f_s = tl.get(fwd["label"], 0)
            b_s = tl.get("NO + N --> N2 + O", 0)
            x = composition(run, RES_SPC)
            c = x["N2"] * x["O"] / (x["NO"] * x["N"])
            Bc += b_s * c
            per.append(f_s / max(b_s * c, 1))
            kden += (x["NO"] * NLOW) * (x["N"] * NLOW) * V * NST * DT / fnum
        r, e = ratio_err(F, B, per)
        r, e = F / Bc, e * (F / Bc) / r
        kb = B / kden
        kb_e = kb / math.sqrt(B)
        capf = XR.mean_p(T) / XR.mean_p(T, cap=False)
        zel[Tn] = (r / keq, e / keq, kb, kb_e, capf)
        num["zeld_ratio_%d" % Tn] = [r / keq, 2 * e / keq, 1.0 / capf, T]
        kd = R.arrhenius(fwd, T) / keq
        num["zeld_kb_%d" % Tn] = [kb, 2 * kb_e, kd, kd * capf, T]
        rep.check(abs(r / keq * capf - 1) <= 3 * e / keq + 0.005,
                  "Zeldovich T=%5.0f K: (f/b)/K_eq = %.3f +/- %.3f (%d/%d "
                  "tallies, realized composition), expected %.4f (reverse "
                  "capped to %.4f of k_f/K_eq)" % (T, r / keq, 2 * e / keq,
                                                   F, B, 1 / capf, capf))
        rep.check(abs(kb / (kd * capf) - 1) <= 3 * kb_e / kd + 0.005,
                  "Zeldovich T=%5.0f K: reverse rate %.4e +/- %.1e, derived "
                  "k_f/K_eq x capped fraction %.4e" % (T, kb, 2 * kb_e,
                                                       kd * capf))
    # legacy published pair
    leg = R.parse_tce(open(os.path.join(HERE, "legacy_air.tce")).read(),
                      keep_reverse=True)
    lf = [x for x in leg if x["label"] == "N2 + O --> NO + N"][0]
    lb = [x for x in leg if x["label"] == "NO + N --> N2 + O"][0]
    Tg = np.linspace(7000.0, 16000.0, 10)
    imp = [R.arrhenius(lf, T) / R.arrhenius(lb, T) / D.keq(fwd, T) for T in Tg]
    der = [R.arrhenius(fwd, T) / D.keq(fwd, T) / R.arrhenius(lb, T)
           for T in Tg]
    num["legacy_implied_keq_over_pf_7_16kK"] = [min(imp), max(imp)]
    num["derived_over_legacy_kb_7_16kK"] = [min(der), max(der)]
    rep.info("legacy pair: implied K_eq / partition-function K_eq = %.2f-%.2f "
             "over 7-16 kK; derived / legacy reverse rate = %.2f-%.2f"
             % (min(imp), max(imp), min(der), max(der)))

    # ---------------- recombination
    rec = {}
    for name, (tce, vib, rmod, T) in REC.items():
        tot, sel = sums(rl, "rec_" + name)
        if not tot:
            continue
        T = mean_T(sel)
        if "mol" in name:
            fl, bl = "N2 + N2 --> N + N + N2", "N + N --> N2 + N2"
        else:
            fl, bl = "N2 + N --> N + N + N", "N + N --> N2 + N"
        F, B = tot.get(fl, 0), tot.get(bl, 0)
        if F == 0 or B == 0:
            rep.check(False, "recombination %s: tallies %d/%d" % (name, F, B))
            continue
        per = [pkg.tallies(r.log).get(fl, 0) /
               max(pkg.tallies(r.log).get(bl, 0), 1) for r in sel]
        # nominal composition: these long runs (in.res) do not write the
        # species counts.  The sampled x_N^2/x_N2 differs from the nominal
        # value by about 2% per seed, below the counting error of the
        # ratio (3-5%); the seed spread in ratio_err includes it.
        r, e = ratio_err(F, B, per)
        nN, nN2 = 0.25 * NHIGH, 0.25 * NHIGH
        meas = r * nN * nN / nN2
        keq = D.keq(D.reaction(fl), T)
        if vib == "smooth":   # classical vibration of N2 in the K_eq
            th = D.sp["N2"]["vibtemp"]
            keq *= (1.0 / (1.0 - math.exp(-th / T))) / (T / th)
        dev, res = meas / keq - 1, e / r
        rec[name] = (dev, res, F, B, T)
        num["rec_%s" % name] = [dev, 2 * res]
        rep.check(abs(dev) <= 3 * res + 0.01,
                  "recombination %-9s T=%5.0f K: K_eq measured/partition "
                  "function = %.3f +/- %.3f (%d/%d tallies)"
                  % (name, T, 1 + dev, 2 * res, F, B))

    # ---------------- full mechanism
    mech = {}
    for T in (15000.0, 8000.0):
        tot, sel = sums(rl, "mech_T%d" % T)
        if not tot:
            continue
        mech[T] = _mechanism(D, sel, sel[0].variables["NRHO"])
        lines = ["%-28s %7s %8s %8s %7s %6s" % ("pair (forward)", "type",
                                                 "f", "b", "ratio", "2sd")]
        nin = nall = 0
        for lab, typ, F, B, ratio, e in mech[T]:
            if ratio is None:
                lines.append("%-28s %7s %8d %8d   (too few tallies)"
                             % (lab[:28], typ, F, B))
                continue
            lines.append("%-28s %7s %8d %8d %7.3f %5.1f%%"
                         % (lab[:28], typ, F, B, ratio, 200 * e))
            if typ in ("D", "E"):
                nall += 1
                nin += abs(ratio - 1) <= 3 * e + 0.01
        open(os.path.join(out, "mechanism_%dK.txt" % T), "w").write(
            "\n".join(lines) + "\n")
        num["mech_%d_within" % T] = [nin, nall]
        rep.info("mechanism at %5d K (measured %.0f K): %d of %d "
                 "dissociation/exchange pairs within 3 sd of detailed "
                 "balance (see mechanism_%dK.txt)"
                 % (T, np.mean([run_T(q) for q in sel]), nin, nall, T))

    _fig7(out, zel, D, fwd, lf, lb)
    _table7(out, rec)
    if mech:
        _figmech(out, mech)
    with open(os.path.join(out, "numbers.json"), "w") as fh:
        json.dump(num, fh, indent=1, default=float)
    return rep.write()


def _mechanism(D, sel, ntot):
    """(label, type, f, b, measured/predicted, error) per forward/reverse
    pair, pooled over the runs sel.  Each run's prediction of f/b uses its
    own composition x (species counts) and time-averaged temperature T:
    exchange f/b = K_eq(T) x_A x_B/(x_C x_D); dissociation (per third
    body) f/b = K_eq(T) x_AB/(x_a x_b ntot), K_eq in number-density
    units, ntot the total number density.  The measured/predicted ratio is
    sum(f)/sum(b * prediction); its error is the counting error."""
    fwd = [r for r in D.rx_all if r["style"] == "A"]
    rev = [r for r in D.rx_all if r["style"] == "B"]
    per = [(pkg.tallies(q.log), composition(q, MECH_SPC), run_T(q))
           for q in sel]
    rows = []
    for f in fwd:
        if f["type"] == "D":
            AB, M = f["reac"]
            a, b = f["prod"][0], f["prod"][1]
            r = [x for x in rev if x["reac"][:2] in ([a, b], [b, a])
                 and x["prod"] == [AB, M]]
            pred = lambda n, T: D.keq(f, T) * n[AB] / (n[a] * n[b] * ntot)
        elif f["type"] in ("E", "I"):
            r = [x for x in rev if sorted(x["reac"]) == sorted(f["prod"])
                 and sorted(x["prod"]) == sorted(f["reac"])]
            A, B = f["reac"]
            C, Dd = f["prod"][0], f["prod"][1]
            pred = lambda n, T: D.keq(f, T) * n[A] * n[B] / (n[C] * n[Dd])
        else:
            continue
        if not r:
            continue
        F = sum(t.get(f["label"], 0) for t, _, _ in per)
        B_ = sum(t.get(r[0]["label"], 0) for t, _, _ in per)
        if F < 5 or B_ < 5:
            rows.append((f["label"], f["type"], F, B_, None, None))
            continue
        W = sum(t.get(r[0]["label"], 0) * pred(n, T) for t, n, T in per)
        ratio = F / W
        rows.append((f["label"], f["type"], F, B_, ratio,
                     ratio * math.sqrt(1.0 / F + 1.0 / B_)))
    return rows


def _setup():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import paper_style as ps
    ps.apply()
    return plt, ps


def _fig7(out, zel, D, fwd, lf, lb):
    plt, ps = _setup()
    fig, axs = plt.subplots(1, 3, figsize=(7.0, 2.6))
    a, b, c = axs
    Tg = np.linspace(6000.0, 17000.0, 60)
    a.plot(Tg / 1e3, [R.arrhenius(fwd, T) / D.keq(fwd, T) for T in Tg],
           color=ps.BLUE, label="derived, $k_f/K_{\\mathrm{eq}}$")
    a.plot(Tg / 1e3, [R.arrhenius(lb, T) for T in Tg], color=ps.VERM,
           ls="--", label="legacy fit")
    pts = [(T, zel[T][2], zel[T][3]) for T in ZELD_T if T in zel]
    if pts:
        x, y, e = zip(*pts)
        a.errorbar(np.array(x) / 1e3, y, yerr=2 * np.array(e), fmt="o",
                   color=ps.BLUE, mfc="white", ms=4.5, capsize=2, lw=0,
                   label="SPARTA")
    a.set_yscale("log")
    a.set_xlabel("$T$ (kK)")
    a.set_ylabel("$k_r$ (m$^3$ s$^{-1}$)")
    a.legend(fontsize=6.5)
    ps.panel(a, "(a)", "se")
    Tb = np.linspace(2000.0, 20000.0, 80)
    b.plot(Tb / 1e3, [R.arrhenius(lf, T) / R.arrhenius(lb, T) / D.keq(fwd, T)
                      for T in Tb], color=ps.VERM)
    b.axhline(1.0, color=ps.GRAY, lw=0.8, ls=":")
    b.set_xlabel("$T$ (kK)")
    b.set_ylabel("$K_{\\mathrm{eq}}^{\\mathrm{fits}}/K_{\\mathrm{eq}}$")
    ps.panel(b, "(b)", "se")
    pts = [(T, zel[T][0], zel[T][1]) for T in ZELD_T if T in zel]
    if pts:
        x, y, e = zip(*pts)
        c.errorbar(np.array(x) / 1e3, y, yerr=2 * np.array(e), fmt="o",
                   color=ps.BLUE, mfc="white", ms=4.5, capsize=2, lw=0)
    c.axhline(1.0, color=ps.GRAY, lw=0.8, ls=":")
    c.set_xlabel("$T$ (kK)")
    c.set_ylabel("$(N_f/N_r)/K_{\\mathrm{eq}}$")
    ps.panel(c, "(c)", "se")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig02_zeldovich.pdf"))
    plt.close(fig)


def _table7(out, rec):
    names = (("atom", "atomic third body, discrete vibration"),
             ("keq", "equilibrium-constant fit applied to the pair"),
             ("mol", "molecular third body (N$_2$)"),
             ("molsmooth", "molecular third body, continuous vibration"),
             ("atom6765", "atomic third body"))
    L = ["configuration & $T$ (K) & deviation & resolution \\\\",
         "\\hline"]
    for k, lab in names:
        if k in rec:
            dev, res, F, B, T = rec[k]
            L.append("%s & %d & %.1f\\%% & %.1f\\%% \\\\"
                     % (lab, T, 100 * abs(dev), 200 * res))
    open(os.path.join(out, "table7_recomb.tex"), "w").write("\n".join(L) + "\n")
    open(os.path.join(out, "table7_recomb.txt"), "w").write(
        "\n".join(l.replace("\\\\", "").replace("&", "|") for l in L) + "\n")


def _figmech(out, mech):
    plt, ps = _setup()
    fig, axs = plt.subplots(len(mech), 1, figsize=(7.0, 2.4 * len(mech)),
                            sharex=True)
    axs = np.atleast_1d(axs)
    for ax, (T, rows) in zip(axs, sorted(mech.items(), reverse=True)):
        labs = [r[0] for r in rows]
        for i, (lab, typ, F, B, ratio, e) in enumerate(rows):
            if ratio is None:
                continue
            col = {"D": ps.BLUE, "E": ps.ORANGE, "I": ps.GREEN}[typ]
            ax.errorbar(i, ratio, yerr=2 * e, fmt="o", color=col, mfc="white",
                        ms=4, capsize=1.5, lw=0)
        ax.axhline(1.0, color=ps.GRAY, lw=0.8, ls=":")
        ax.set_ylabel("measured/predicted")
        ax.set_yscale("log")
        ax.text(0.01, 0.95, "%d K" % T, transform=ax.transAxes, va="top",
                fontsize=8)
    axs[-1].set_xticks(range(len(labs)))
    axs[-1].set_xticklabels([l.replace(" --> ", r"$\to$") for l in labs],
                            rotation=90, fontsize=5.5)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig_mechanism.pdf"))
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(analyze(sys.argv[1] if len(sys.argv) > 1 else "quick"))
