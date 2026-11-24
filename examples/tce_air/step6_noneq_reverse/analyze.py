"""Step 6 analysis: reverse rates without a cell temperature.

Produces, in results/step6_noneq_reverse/:
  table_cellT.tex/.txt  rate relative to its value at
                          T_int = T_tr = 8000 K, measured and predicted, for
                          the reverse of each configuration (and the forward
                          of data/air.tce), and the value of a model
                          evaluated at the translational cell temperature
  fig_cellT.pdf           the same as a figure: (a) all internal modes at
                          T_int; (b) one mode at a time, and three-body
                          recombination against T_el
  numbers.json, results.txt

References (tce_ref.py): ExchangeReverse, the derived reverse probability
of Eq. (exch) with its absolute calibration and the cap at one, averaged
over the collision ensemble of the reverse pair with each mode at its own
temperature; the same construction for the forward microcanonical
probability (forward=True); Recombination3, the three-body average of
Eq. (pr3) with the electronic ladders at T_el.
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
from runs import (runs, STEP, TTR, TINT, MODES, MODE_T, TEL,  # noqa: E402
                  PAIRS)

LAB = {"rel": ("N2 + O --> NO + N", "NO + N --> N2 + O"),
       "bar": ("NO + N --> N2 + O", "N2 + O --> NO + N")}


def counts(rl, prefix, label):
    sel = [r for r in rl if r.tag.startswith(prefix + "_s") and r.done()]
    per = [pkg.tallies(r.log).get(label, 0) for r in sel]
    return sum(per), per


def ratio(rl, pre, pre0, label):
    """tally ratio of two run sets (same particles and steps) with its
    error: Poisson, or the spread over seeds if larger"""
    n, per = counts(rl, pre, label)
    n0, per0 = counts(rl, pre0, label)
    if n == 0 or n0 == 0:
        return None, None
    r = n / n0
    e = r * math.sqrt(1.0 / n + 1.0 / n0)
    if len(per) > 1 and len(per) == len(per0):
        rs = [a / b for a, b in zip(per, per0) if b > 0]
        if len(rs) > 1:
            e = max(e, np.std(rs, ddof=1) / math.sqrt(len(rs)))
    return r, e


def check_frozen_files(rep):
    def rows(path, skip):
        out = []
        for l in open(path):
            w = l.split()
            if l.strip() and not l.lstrip().startswith("#") and len(w) > 2:
                out.append([x for i, x in enumerate(w) if i not in skip(w)])
        return out
    sp_skip = lambda w: {4, 6}
    el_skip = lambda w: {2 + 5 * k + 1 for k in range(int(w[1]))}
    for f, ref, sk in (("air_frozen.species", "air.species", sp_skip),
                       ("air_frozen.elec", "air.elec", el_skip)):
        a = rows(os.path.join(HERE, f), sk)
        b = rows(os.path.join(R.DATA, ref), sk)
        same = len(a) == len(b) and all(
            len(x) == len(y) and all(float(p) == float(q) if
                                     _isnum(p) else p == q
                                     for p, q in zip(x, y))
            for x, y in zip(a, b))
        rep.check(same, "%s equals data/%s except the relaxation "
                  "probabilities" % (f, ref))


def _isnum(x):
    try:
        float(x)
        return True
    except ValueError:
        return False


def analyze(size="quick"):
    rep = pkg.Report(STEP + " (" + size + ")")
    out = os.path.join(pkg.results_dir(), STEP)
    os.makedirs(out, exist_ok=True)
    rl = runs(size)
    check_frozen_files(rep)
    D = R.Data()
    num, table = {}, []
    refs = {}
    for pair, tce in PAIRS.items():
        Dp = R.Data(tce=os.path.join(HERE, tce))
        fwd = Dp.reaction(LAB[pair][0])
        refs[pair] = (R.ExchangeReverse(Dp, fwd),
                      R.ExchangeReverse(Dp, fwd, forward=True))

    # ---- all internal modes at T_int
    for pair in PAIRS:
        rev, fw = refs[pair]
        base_r, base_f = rev.mean_p(TTR), fw.mean_p(TTR)
        for Ti in TINT:
            pre, pre0 = "tint_%s_TI%d" % (pair, Ti), "tint_%s_TI%d" % (
                pair, TTR)
            pr = rev.mean_p(TTR, Trot=Ti, Tvib=Ti, Telec=Ti) / base_r
            pf = fw.mean_p(TTR, Trot=Ti, Tvib=Ti, Telec=Ti) / base_f
            mr, er = ratio(rl, pre, pre0, LAB[pair][1])
            mf, ef = ratio(rl, pre, pre0, LAB[pair][0])
            if mr is None:
                continue
            table.append((pair, Ti, mr, er, pr, mf, ef, pf))
            num["tint_%s_%d" % (pair, Ti)] = dict(
                reverse=[mr, 2 * er, pr], forward=[mf, 2 * ef, pf])
            if Ti == TTR:
                continue
            rep.check(abs(mr / pr - 1) <= 3 * er / pr + 0.005,
                      "%s T_int=%5d K: reverse %-18s x%.3f +/- %.3f, "
                      "predicted x%.3f" % (pair, Ti, LAB[pair][1], mr,
                                           2 * er, pr),
                      meas=mr, ref=pr, sd=er, tol=0.005 * pr,
                      label="%s, T_int %d K, reverse" % (pair, Ti))
            rep.check(abs(mf / pf - 1) <= 3 * ef / pf + 0.005,
                      "%s T_int=%5d K: forward %-18s x%.3f +/- %.3f, "
                      "predicted x%.3f" % (pair, Ti, LAB[pair][0], mf,
                                           2 * ef, pf),
                      meas=mf, ref=pf, sd=ef, tol=0.005 * pf,
                      label="%s, T_int %d K, forward" % (pair, Ti))

    # ---- one mode at a time (rel pair)
    rev, fw = refs["rel"]
    base_r = rev.mean_p(TTR)
    for mode, var in MODES.items():
        for lab, Tm in MODE_T.items():
            kw = {"Trot": TTR, "Tvib": TTR, "Telec": TTR}
            kw[{"TROT": "Trot", "TVIB": "Tvib", "TELEC": "Telec"}[var]] = Tm
            pr = rev.mean_p(TTR, **kw) / base_r
            mr, er = ratio(rl, "mode_%s_%s" % (mode, lab),
                           "tint_rel_TI%d" % TTR, LAB["rel"][1])
            if mr is None:
                continue
            num["mode_%s_%s" % (mode, lab)] = [mr, 2 * er, pr]
            rep.check(abs(mr / pr - 1) <= 3 * er / pr + 0.005,
                      "%-4s at %5d K, others at %d K: reverse x%.3f +/- "
                      "%.3f, predicted x%.3f" % (mode, Tm, TTR, mr, 2 * er,
                                                 pr),
                      meas=mr, ref=pr, sd=er, tol=0.005 * pr)

    # ---- recombination with the electronic ladders at T_el
    rc = R.Recombination3(D, D.reaction("N2 + N --> N + N + N"))
    for Te in TEL:
        if Te == TTR:
            continue
        mr, er = ratio(rl, "recel_TE%d" % Te, "recel_TE%d" % TTR,
                       "N + N --> N2 + N")
        if mr is None:
            continue
        pr = rc.ratio(TTR, Te)
        num["recel_%d" % Te] = [mr, 2 * er, pr]
        rep.check(abs(mr / pr - 1) <= 3 * er / pr + 0.005,
                  "recombination N + N + N, T_el=%5d K: x%.3f +/- %.3f, "
                  "three-body quadrature x%.4f" % (Te, mr, 2 * er, pr),
                  meas=mr, ref=pr, sd=er, tol=0.005 * pr)
    _table(out, table)
    _fig(out, table, num, refs)
    with open(os.path.join(out, "numbers.json"), "w") as fh:
        json.dump(num, fh, indent=1, default=float)
    return rep.write()


def _fig(out, table, num, refs):
    """fig_cellT.pdf: (a) reverse and forward rate against T_int, all
    internal modes at T_int and translation at T_tr; (b) one mode at a
    time, and three-body recombination against T_el.  Symbols: SPARTA
    with 2 sd; lines and bars: tce_ref.  A cell-temperature model gives 1
    everywhere (dashed)."""
    if not table:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import paper_style as ps
    ps.apply()
    fig, (ax, bx) = plt.subplots(1, 2, figsize=ps.DOUBLE,
                                 gridspec_kw=dict(width_ratios=(1.15, 1)))
    Tg = np.linspace(3000.0, 22000.0, 20)
    series = (("rel", 1, ps.BLUE, "o", "derived NO + N \u2192 N$_2$ + O"),
              ("rel", 0, ps.ORANGE, "s", "fitted N$_2$ + O \u2192 NO + N"),
              ("bar", 1, ps.GREEN, "^", "derived N$_2$ + O \u2192 NO + N"),
              ("bar", 0, ps.PURPLE, "v", "fitted NO + N \u2192 N$_2$ + O"))
    for pair, isrev, col, mk, lab in series:
        rev, fw = refs[pair]
        ref = rev if isrev else fw
        base = ref.mean_p(TTR)
        ax.plot(Tg, [ref.mean_p(TTR, Trot=t, Tvib=t, Telec=t) / base
                     for t in Tg], color=col, lw=1.2)
        pts = [(Ti, mr, er) if isrev else (Ti, mf, ef)
               for p_, Ti, mr, er, pr, mf, ef, pf in table if p_ == pair]
        ax.errorbar([x[0] for x in pts], [x[1] for x in pts],
                    yerr=[2 * x[2] for x in pts], fmt=mk, color=col,
                    mfc="white", capsize=1.5, lw=0, elinewidth=0.9,
                    label=lab)
    ax.axhline(1.0, color=ps.GRAY, lw=0.9, ls="--",
               label="cell-temperature $k_f/K_{eq}(T_{tr})$")
    ax.axvline(TTR, color=ps.GRAY, lw=0.5, ls=":")
    ax.set_yscale("log")
    ax.set_xlabel("$T_{int}$ (K), all internal modes; $T_{tr}$ = %d K"
                  % TTR)
    ax.set_ylabel("rate / rate at $T_{int} = T_{tr}$")
    ax.legend(loc="upper left", fontsize=6.5)
    ps.panel(ax, "(a)", loc="se")
    cats, meas, err, pred = [], [], [], []
    for mode in MODES:
        for lab in MODE_T:
            k = "mode_%s_%s" % (mode, lab)
            if k in num:
                cats.append("%s %s" % (mode, lab))
                meas.append(num[k][0]); err.append(num[k][1])
                pred.append(num[k][2])
    for Te in TEL:
        k = "recel_%d" % Te
        if k in num:
            cats.append("recomb.\n$T_{el}$ %d" % Te)
            meas.append(num[k][0]); err.append(num[k][1])
            pred.append(num[k][2])
    x = np.arange(len(cats))
    bx.bar(x, pred, width=0.6, color="0.85", label="tce_ref")
    bx.errorbar(x, meas, yerr=err, fmt="o", color=ps.BLUE, mfc="white",
                capsize=1.5, lw=0, elinewidth=0.9, label="SPARTA (2 sd)")
    bx.axhline(1.0, color=ps.GRAY, lw=0.9, ls="--")
    bx.set_xticks(x)
    bx.set_xticklabels(cats, rotation=60, fontsize=6.5)
    bx.minorticks_off()
    bx.set_ylabel("reverse rate / rate at $T_{tr}$")
    bx.legend(loc="upper left", fontsize=6.5)
    ps.panel(bx, "(b)", loc="ne")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig_cellT.pdf"))
    plt.close(fig)


def _table(out, table):
    L = ["configuration & $T_{\\mathrm{int}}$ (K) & reverse measured & "
         "microcanonical & forward measured & microcanonical & "
         "translational cell temperature \\\\", "\\hline"]
    T = ["%-4s %8s %10s %8s %10s %8s %6s" % ("pair", "T_int", "rev meas",
                                             "rev pred", "fwd meas",
                                             "fwd pred", "cell")]
    for pair, Ti, mr, er, pr, mf, ef, pf in table:
        L.append("%s & %d & %.3f & %.3f & %.3f & %.3f & 1.000 \\\\"
                 % ({"rel": "air.tce", "bar": "reversed"}[pair], Ti, mr, pr,
                    mf, pf))
        T.append("%-4s %8d %10.3f %8.3f %10.3f %8.3f %6.3f"
                 % (pair, Ti, mr, pr, mf, pf, 1.0))
    open(os.path.join(out, "table_cellT.tex"), "w").write("\n".join(L) + "\n")
    open(os.path.join(out, "table_cellT.txt"), "w").write("\n".join(T) + "\n")


if __name__ == "__main__":
    sys.exit(analyze(sys.argv[1] if len(sys.argv) > 1 else "quick"))
