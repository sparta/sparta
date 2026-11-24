"""Step 4 analysis: closed reacting box with the Zeldovich pair.

Produces, in results/step4_box/:
  fig_closedbox.pdf  NO fraction of the NO + N2 pool against time,
                       (a) thermostatted at 15 kK, (b) adiabatic, from both
                       sides, with the equilibrium values
  numbers.json, results.txt (with the mode temperatures of each species)

References (tce_ref.py, from the data files only):
  isothermal: K_eq(T) = n_NO n_N/(n_N2 n_O); with n_N2 = n_O and
              n_NO = n_N the NO fraction is sqrt(K)/(1 + sqrt(K))
  adiabatic:  the composition and temperature at which K_eq holds and the
              total energy (thermal energy of all modes plus formation
              energy of data/air.tce) equals that of the initial sample
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
from runs import runs, STEP, T0, STARTS  # noqa: E402

LABEL = "N2 + O --> NO + N"


def series(run):
    """(time in us, NO/(NO+N2), translational T, counts at step 0)"""
    head, rows = pkg.stats(run.log)
    a = np.array(rows)
    col = {h: i for i, h in enumerate(head)}
    n2, no = a[:, col["c_spc[1]"]], a[:, col["c_spc[3]"]]
    counts0 = {s: a[0, col["c_spc[%d]" % (i + 1)]]
               for i, s in enumerate(("N2", "O", "NO", "N"))}
    return a[:, col["Step"]] * 1e-3, no / (no + n2), a[:, col["c_temp"]], \
        counts0


MODES = (("T_tr", "c_temp"), ("T_rot N2", "c_trotave[1]"),
         ("T_rot NO", "c_trotave[2]"), ("T_vib N2", "c_tvibave[1]"),
         ("T_vib NO", "c_tvibave[2]"), ("T_el N2", "c_telave[1]"),
         ("T_el O", "c_telave[2]"), ("T_el NO", "c_telave[3]"),
         ("T_el N", "c_telave[4]"))


def mode_temps(sel):
    """mode temperatures averaged over the second half and over seeds"""
    acc = {}
    for r in sel:
        head, rows = pkg.stats(r.log)
        a = np.array(rows)
        col = {h: i for i, h in enumerate(head)}
        st = a[:, col["Step"]]
        keep = st >= 0.5 * st[-1]
        for lab, c in MODES:
            if c in col:
                acc.setdefault(lab, []).append(a[keep, col[c]].mean())
    return [(lab, float(np.mean(acc[lab]))) for lab, c in MODES if lab in acc]


def tail_mean(t, x, nblock=5):
    """mean over the second half of the run and its error from block means"""
    sel = x[t >= 0.5 * t[-1]]
    blocks = np.array_split(sel, nblock)
    m = np.array([b.mean() for b in blocks])
    return sel.mean(), m.std(ddof=1) / math.sqrt(nblock)


def adiabatic_eq(D, fwd, counts0, T0):
    """equilibrium temperature and NO fraction of the closed adiabatic box
    with the initial composition counts0 (particle counts) at T0"""
    sp = ("N2", "O", "NO", "N")
    E0 = sum(counts0[s] * (D.e_mean(s, T0) + R.EFORM[s]) for s in sp)

    def comp(T):
        # N2 = a - x, O = b - x, NO = c + x, N = d + x with K_eq(T)
        a, b, c, d = (counts0[s] for s in sp)
        K = D.keq(fwd, T)
        f = lambda x: (c + x) * (d + x) - K * (a - x) * (b - x)
        lo, hi = -min(c, d), min(a, b)
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if f(mid) > 0:
                hi = mid
            else:
                lo = mid
        x = 0.5 * (lo + hi)
        return {"N2": a - x, "O": b - x, "NO": c + x, "N": d + x}

    def excess(T):
        n = comp(T)
        return sum(n[s] * (D.e_mean(s, T) + R.EFORM[s]) for s in sp) - E0

    lo, hi = 2000.0, 60000.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if excess(mid) > 0:
            hi = mid
        else:
            lo = mid
    T = 0.5 * (lo + hi)
    n = comp(T)
    return T, n["NO"] / (n["NO"] + n["N2"])


def analyze(size="quick"):
    rep = pkg.Report(STEP + " (" + size + ")")
    out = os.path.join(pkg.results_dir(), STEP)
    os.makedirs(out, exist_ok=True)
    D = R.Data()
    fwd = D.reaction(LABEL)
    K = D.keq(fwd, T0)
    x_iso = math.sqrt(K) / (1.0 + math.sqrt(K))
    num = {"keq_15kK": K, "x_iso": x_iso}
    rep.info("isothermal equilibrium at %.0f K: K_eq = %.4f, NO fraction "
             "%.4f" % (T0, K, x_iso))
    curves = {}
    for mode in ("iso", "adi"):
        for side in STARTS:
            sel = [r for r in runs(size)
                   if r.tag.startswith("%s_%s_s" % (mode, side)) and r.done()]
            if not sel:
                continue
            means, errs, refs, Teq, Tm = [], [], [], [], []
            for r in sel:
                t, x, Ttr, c0 = series(r)
                curves.setdefault((mode, side), []).append((t, x))
                m, e = tail_mean(t, x)
                means.append(m)
                errs.append(e)
                Tm.append(tail_mean(t, Ttr)[0])
                if mode == "adi":
                    Te, xe = adiabatic_eq(D, fwd, c0, T0)
                    refs.append(xe)
                    Teq.append(Te)
                else:
                    refs.append(x_iso)
            mt = mode_temps(sel)
            rep.info("%s from the %s side, mean over the second half: " %
                     ("thermostatted" if mode == "iso" else "adiabatic",
                      side) + ", ".join("%s %.0f K" % kv for kv in mt))
            num["%s_%s_modeT" % (mode, side)] = dict(mt)
            m = float(np.mean(means))
            e = math.sqrt(sum(v * v for v in errs)) / len(errs)
            if len(means) > 1:
                e = max(e, np.std(means, ddof=1) / math.sqrt(len(means)))
            ref = float(np.mean(refs))
            num["%s_%s" % (mode, side)] = [m, 2 * e, ref]
            x0 = curves[(mode, side)][0][1][0]
            if mode == "iso":
                rep.check(abs(m - ref) <= 3 * e + 0.003,
                          "thermostatted, from %.3f: NO fraction %.4f +/- "
                          "%.4f (2 sd), equilibrium %.4f"
                          % (x0, m, 2 * e, ref),
                          meas=m, ref=ref, sd=e, tol=0.003)
            else:
                num["adi_%s_T" % side] = [float(np.mean(Tm)),
                                          float(np.mean(Teq))]
                rep.check(abs(m - ref) <= 3 * e + 0.005,
                          "adiabatic, from %.3f: NO fraction %.4f +/- %.4f "
                          "(2 sd), energy-conserving equilibrium %.4f "
                          "(T %.0f K measured, %.0f K predicted)"
                          % (x0, m, 2 * e, ref, np.mean(Tm), np.mean(Teq)),
                          meas=m, ref=ref, sd=e, tol=0.005)
    _fig(out, curves, num)
    with open(os.path.join(out, "numbers.json"), "w") as fh:
        json.dump(num, fh, indent=1, default=float)
    return rep.write()


def _fig(out, curves, num):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, os.path.dirname(HERE))
    import paper_style as ps
    ps.apply()
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.6), sharey=True)
    for ax, mode, title in zip(axs, ("iso", "adi"),
                               ("thermostatted at 15 kK", "adiabatic")):
        for side, col in (("fwd", ps.BLUE), ("bwd", ps.VERM)):
            for i, (t, x) in enumerate(curves.get((mode, side), [])):
                ax.plot(t, x, color=col, lw=0.9, alpha=1.0 if i == 0 else 0.4,
                        label=("from N$_2$ + O" if side == "fwd"
                               else "from NO + N") if i == 0 else None)
            key = "%s_%s" % (mode, side)
            if key in num:
                ax.axhline(num[key][2], color=col if mode == "adi" else
                           ps.GRAY, ls="--", lw=0.8)
        ax.set_xlabel("time ($\\mu$s)")
        ax.text(0.97, 0.5, title, transform=ax.transAxes, ha="right",
                va="center", fontsize=7)
    axs[0].set_ylabel("NO/(NO + N$_2$)")
    axs[0].legend(fontsize=6.5, loc="upper right")
    ps.panel(axs[0], "(a)", "nw")
    ps.panel(axs[1], "(b)", "nw")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig_closedbox.pdf"))
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(analyze(sys.argv[1] if len(sys.argv) > 1 else "quick"))
