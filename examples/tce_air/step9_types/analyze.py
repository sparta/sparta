"""Step 9 analysis: reaction types beyond dissociation.

Produces, in results/step9_types/:
  fig_types.pdf          top: forward equilibrium rate relative to the input
                         Arrhenius rate, per reaction type, for the
                         microcanonical factor and the earlier treatments;
                         bottom: derived reverse rate relative to k_f/K_eq
  table_types.tex/.txt   the same numbers
  numbers.json, results.txt, checks.json, summary.pdf

References (tce_ref.py, from the data files only):
  micro forward   the capped fraction <min(P,1)>/<P> of the forward
                  microcanonical probability (ExchangeReverse, forward=True),
                  1 where nothing is capped
  incl            R_inclusion: electronic energy added to the collision
                  energy, Eq. (5) of the paper, vibration left out
  inst            R_instdof: Bird's instantaneous vibrational DOF, no
                  electronic states
  reverse         k_f/K_eq times the capped fraction of the derived reverse
                  probability (ExchangeReverse), 1 where nothing is capped;
                  dissociative recombination is capped below about 15 kK
                  (the low-temperature limit of the paper)
Every reference is evaluated at the run's measured mean temperature, and
the rates use the composition the sample actually has.

The error of each mean is the counting error or the spread over seeds,
whichever is larger.  For the ionization pair the spread wins by an order
of magnitude (1-1.5% against 0.1%): as in step7_cap, a frozen reservoir
keeps its initial sample of spin states (non-reacting collisions conserve
spin and no reaction is performed), for example the NO+ triplets, about 7%
of NO+ at 20 kK, which carry their own share of the recombination.
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
from runs import runs, STEP, PAIRS, TYPES_T, treatments  # noqa: E402

V = (1.0e-4) ** 3
DT = 1.0e-9
TOL = 0.01


def stats_array(run):
    head, rows = pkg.stats(run.log)
    a = np.array(rows)
    return {h: a[:, i] for i, h in enumerate(head)}


def measure(run, label, species):
    """rate coefficient of reaction label from its tally, the realized
    composition and the run time; with the relative counting error and the
    mean temperature"""
    st = stats_array(run)
    Tm = float(np.mean(st["c_temp"][1:]))
    cnt = [st["c_spc[%d]" % (i + 1)][0] for i in range(4)]
    x = dict(zip(_species(run), [c / sum(cnt) for c in cnt]))
    n = run.variables["NRHO"]
    fnum = n * V / run.variables["NP"]
    t = run.variables["NSTEP"] * DT
    F = pkg.tallies(run.log).get(label, 0)
    if F == 0:
        return None
    a, b = species
    k = F * fnum / (x[a] * n * x[b] * n * V * t)
    return k, 1.0 / math.sqrt(F), Tm, F


def _species(run):
    return [run.variables["S%d" % (i + 1)] for i in range(4)]


def ladder(D, s):
    return D.elec.get(s) or [(0.0, 1.0)]


def mean_err(vals):
    """mean of per-seed ratios; error from counting or the seed spread,
    whichever is larger"""
    y = [v[0] for v in vals]
    m = float(np.mean(y))
    e = math.sqrt(sum((v[0] * v[1]) ** 2 for v in vals)) / len(y)
    if len(y) > 1:
        e = max(e, float(np.std(y, ddof=1)) / math.sqrt(len(y)))
    return m, e


def analyze(size="quick"):
    rep = pkg.Report(STEP + " (" + size + ")")
    out = os.path.join(pkg.results_dir(), STEP)
    os.makedirs(out, exist_ok=True)
    rl = runs(size)
    D = R.Data()
    num, rows = {}, []
    refs = {}
    for pair, (sp, x, flab, blab, lab) in PAIRS.items():
        fwd = D.reaction(flab)
        A, B = fwd["reac"]
        z0, om = D.zrot(A, B), D.pair(A, B)[1]
        th = (D.sp[A]["vibtemp"], D.sp[B]["vibtemp"])
        fw = R.ExchangeReverse(D, fwd, forward=True)
        rv = R.ExchangeReverse(D, fwd)
        refs[pair] = dict(
            micro=lambda T, fw=fw: fw.mean_p(T) / fw.mean_p(T, cap=False),
            incl=lambda T, A=A, B=B, z0=z0, om=om, f=fwd: R.R_inclusion(
                T, ladder(D, A), ladder(D, B), z0, om, f["eta"], f["Ea"]),
            inst=lambda T, th=th, z0=z0, om=om, f=fwd: R.R_instdof(
                T, th[0], th[1], z0, om, f["eta"], f["Ea"]),
            rev=lambda T, rv=rv: rv.mean_p(T) / rv.mean_p(T, cap=False))
        rep.info("%s: %s (Ea %.3f eV, eta %+.2f); reverse %s"
                 % (lab, flab, fwd["Ea"] / R.EV, fwd["eta"], blab))
        for tr in treatments(pair):
            for T in TYPES_T:
                sel = [r for r in rl if r.tag.startswith(
                    "%s_%s_T%d_s" % (pair, tr, T)) and r.done()]
                if not sel:
                    continue
                f_vals, b_vals, Tms, nF, nB = [], [], [], 0, 0
                for r in sel:
                    mf = measure(r, flab, fwd["reac"])
                    if mf is None:
                        continue
                    k, e, Tm, F = mf
                    Tms.append(Tm)
                    nF += F
                    f_vals.append((k / R.arrhenius(fwd, Tm), e))
                    if tr == "micro":
                        mb = measure(r, blab, fwd["prod"])
                        if mb is not None:
                            kb, eb, _, Bn = mb
                            kdb = R.arrhenius(fwd, Tm) / D.keq(fwd, Tm)
                            b_vals.append((kb / kdb, eb))
                            nB += Bn
                if not f_vals:
                    rep.check(False, "%s %s T=%d K: no forward tallies"
                              % (lab, tr, T))
                    continue
                Tm = float(np.mean(Tms))
                m, e = mean_err(f_vals)
                ref = refs[pair][tr](Tm)
                key = "%s_%s_%d" % (pair, tr, T)
                num[key] = dict(forward=[m, 2 * e, ref], T=Tm, tallies=nF)
                rep.check(abs(m - ref) <= 3 * e + TOL * ref,
                          "%s, %s, T=%5.0f K: k/k_Arr %.3f +/- %.3f (2 sd, "
                          "%d tallies), reference %.3f"
                          % (lab, tr, Tm, m, 2 * e, nF, ref),
                          meas=m, ref=ref, sd=e, tol=TOL * ref,
                          label="%s, %s, %.0f K" % (pair, tr, T))
                row = [pair, tr, Tm, m, e, ref, None, None, None]
                if b_vals:
                    mb, eb = mean_err(b_vals)
                    refb = refs[pair]["rev"](Tm)
                    num[key]["reverse"] = [mb, 2 * eb, refb]
                    num[key]["reverse_tallies"] = nB
                    rep.check(abs(mb - refb) <= 3 * eb + TOL * refb,
                              "%s reverse, T=%5.0f K: k_b/(k_f/K_eq) %.3f "
                              "+/- %.3f (2 sd, %d tallies), capped fraction "
                              "%.3f" % (lab, Tm, mb, 2 * eb, nB, refb),
                              meas=mb, ref=refb, sd=eb, tol=TOL * refb,
                              label="%s, reverse, %.0f K" % (pair, T))
                    row[6:] = [mb, eb, refb]
                rows.append(row)
    _table(out, rows)
    if rows:
        _fig(out, rows, refs)
    with open(os.path.join(out, "numbers.json"), "w") as fh:
        json.dump(num, fh, indent=1, default=float)
    return rep.write()


def _table(out, rows):
    L = ["reaction & treatment & $T$ (K) & $k/k_{\\mathrm{Arr}}$ & "
         "reference & $k_b/(k_f/K_{eq})$ & reference \\\\", "\\hline"]
    T = ["%-5s %-6s %7s %16s %8s %16s %8s" % ("pair", "treat", "T", "k/k_Arr",
                                              "ref", "kb/(kf/Keq)", "ref")]
    for pair, tr, Tm, m, e, ref, mb, eb, refb in rows:
        rv = ("%.3f $\\pm$ %.3f & %.3f" % (mb, 2 * eb, refb)
              if mb is not None else "& ")
        L.append("%s & %s & %.0f & %.3f $\\pm$ %.3f & %.3f & %s \\\\"
                 % (PAIRS[pair][2].replace("-->", "$\\to$"), tr, Tm, m,
                    2 * e, ref, rv))
        T.append("%-5s %-6s %7.0f %9.3f +/-%5.3f %8.3f %s"
                 % (pair, tr, Tm, m, 2 * e, ref,
                    "%9.3f +/-%5.3f %8.3f" % (mb, 2 * eb, refb)
                    if mb is not None else ""))
    open(os.path.join(out, "table_types.tex"), "w").write("\n".join(L) + "\n")
    open(os.path.join(out, "table_types.txt"), "w").write("\n".join(T) + "\n")


def _fig(out, rows, refs):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import paper_style as ps
    ps.apply()
    pairs = list(PAIRS)
    fig, axs = plt.subplots(2, len(pairs), figsize=(7.0, 4.6), sharex=True)
    Tg = np.linspace(min(TYPES_T) * 0.92, max(TYPES_T) * 1.04, 12)
    style = {"micro": (ps.BLUE, "o", "microcanonical (this work)"),
             "incl": (ps.ORANGE, "s", "electronic energy added"),
             "inst": (ps.GREEN, "^", "instantaneous vibrational DOF")}
    for j, pair in enumerate(pairs):
        ax, bx = axs[0, j], axs[1, j]
        trs = sorted({r[1] for r in rows if r[0] == pair},
                     key=list(style).index)
        for tr in trs:
            col, mk, lab = style[tr]
            ax.plot(Tg / 1e3, [refs[pair][tr](t) for t in Tg], color=col,
                    lw=1.1)
            pts = [r for r in rows if r[0] == pair and r[1] == tr]
            ax.errorbar([p[2] / 1e3 for p in pts], [p[3] for p in pts],
                        yerr=[2 * p[4] for p in pts], fmt=mk, color=col,
                        mfc="white", capsize=1.5, lw=0, elinewidth=0.9,
                        label=lab)
        ax.axhline(1.0, color=ps.GRAY, lw=0.7, ls=":")
        flab = PAIRS[pair][2].replace(" --> ", " → ")
        ax.set_title("%s\n%s" % (PAIRS[pair][4], flab), fontsize=7.5)
        if j == 0:
            ax.set_ylabel("$k/k_{\\mathrm{Arr}}$ (forward)")
        bx.plot(Tg / 1e3, [refs[pair]["rev"](t) for t in Tg], color=ps.BLUE,
                lw=1.1, label="capped fraction (tce_ref)")
        pts = [r for r in rows if r[0] == pair and r[6] is not None]
        bx.errorbar([p[2] / 1e3 for p in pts], [p[6] for p in pts],
                    yerr=[2 * p[7] for p in pts], fmt="o", color=ps.BLUE,
                    mfc="white", capsize=1.5, lw=0, elinewidth=0.9,
                    label="SPARTA (2 sd)")
        bx.axhline(1.0, color=ps.GRAY, lw=0.7, ls=":")
        blab = PAIRS[pair][3].replace(" --> ", " → ")
        bx.text(0.04, 0.06, "reverse " + blab, transform=bx.transAxes,
                fontsize=6.5)
        bx.set_xlabel("$T$ (kK)")
        if j == 0:
            bx.set_ylabel("$k_b/(k_f/K_{eq})$ (derived reverse)")
        ps.panel(ax, "(%s)" % "abc"[j], loc="ne")
        ps.panel(bx, "(%s)" % "def"[j], loc="ne")
    axs[0, 0].legend(loc="upper left", fontsize=6)
    axs[1, 0].legend(loc="center left", fontsize=6)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig_types.pdf"))
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(analyze(sys.argv[1] if len(sys.argv) > 1 else "quick"))
