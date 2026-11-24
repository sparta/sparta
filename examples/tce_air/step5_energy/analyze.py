"""Step 5 analysis: energy conservation of the reaction cycles.

Produces, in results/step5_energy/:
  table_energy.tex/.txt  drift of total energy per 100 steps
                          between steps 100 and 300, relative to the initial
                          energy, for the combinations of reaction energies
                          and species masses
  fig_energy.pdf          relative change of total energy against step
  numbers.json, results.txt

Checks: with closed reaction energies and conserved masses the total
energy is constant (|change| < 1e-6 at every output, the paper's
statement); each defect alone produces a drift that is resolved (> 1e-5
per 100 steps).  The remaining change, about 2e-7, is not a defect of
the code: the reaction energies of data/air.tce are written with six
digits and the formation energies of the energy sum with six, and the
two differ in the last digit (2e-24 J per dissociation of N2, for
example), which over the few hundred dissociations of the initial
transient adds up to 1e-7 of the total.
Also checked: the reaction energies of data/air.tce are differences of the
formation-energy table that the energy sum uses (tce_ref.EFORM).
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)
import pkg                 # noqa: E402
import tce_ref as R        # noqa: E402
sys.modules.pop("runs", None)  # this step's runs.py, not another step's
from runs import runs, STEP, CASES  # noqa: E402

ROWS = (("unclosed_unconserved", "not closed", "not conserved"),
        ("unclosed_conserved", "not closed", "conserved"),
        ("closed_unconserved", "closed", "not conserved"),
        ("closed_norecomb", "closed, recombination off", "not conserved"),
        ("closed_conserved", "closed", "conserved"),
        ("full11", "closed (11 species, ions, electrons)", "conserved"))


def energy(run):
    head, rows = pkg.stats(run.log)
    a = np.array(rows)
    col = {h: i for i, h in enumerate(head)}
    step, e = a[:, col["Step"]], a[:, col["v_etot"]]
    return step, e / e[0] - 1.0


def drift(step, de):
    i100, i300 = list(step).index(100), list(step).index(300)
    return (de[i300] - de[i100]) / 2.0


def analyze(size="quick"):
    rep = pkg.Report(STEP + " (" + size + ")")
    out = os.path.join(pkg.results_dir(), STEP)
    os.makedirs(out, exist_ok=True)
    D = R.Data()
    res = max(abs(R.closure_residual(r)) for r in D.rx_all
              if r["style"] == "A")
    rep.check(res < 1e-4, "data/air.tce reaction energies are differences "
              "of the formation-energy table (largest relative residual "
              "%.1e, the rounding of the tabulated digits)" % res)
    num, curves, tab = {}, {}, []
    for name, lab_e, lab_m in ROWS:
        sel = [r for r in runs(size) if r.tag.startswith(name + "_s")
               and r.done()]
        if not sel:
            continue
        ds, mx = [], []
        for r in sel:
            st, de = energy(r)
            curves.setdefault(name, []).append((st, de))
            ds.append(drift(st, de))
            mx.append(float(np.max(np.abs(de))))
        d = float(np.mean(ds))
        num[name] = {"drift_per_100": ds, "max_abs_change": mx}
        tab.append((lab_e, lab_m, d, max(mx)))
        if name in ("closed_conserved", "full11"):
            rep.check(max(mx) < 1e-6,
                      "%-22s total energy constant: largest relative change "
                      "%.1e over the run (drift per 100 steps %+.1e)"
                      % (name, max(mx), d))
        elif name == "closed_norecomb":
            rep.check(abs(d) < 1e-6,
                      "%-22s no reaction cycle: drift per 100 steps %+.1e "
                      "between steps 100 and 300" % (name, d))
        else:
            rep.check(abs(d) > 1e-5,
                      "%-22s drift per 100 steps %+.2e (seeds: %s)"
                      % (name, d, ", ".join("%+.2e" % x for x in ds)))
    if all(k in num for k in ("unclosed_unconserved", "unclosed_conserved",
                              "closed_unconserved")):
        a = np.mean(num["unclosed_unconserved"]["drift_per_100"])
        b = np.mean(num["unclosed_conserved"]["drift_per_100"])
        c = np.mean(num["closed_unconserved"]["drift_per_100"])
        rep.info("additivity: both defects %+.2e, sum of the single defects "
                 "%+.2e" % (a, b + c))
    _table(out, tab)
    _fig(out, curves)
    with open(os.path.join(out, "numbers.json"), "w") as fh:
        json.dump(num, fh, indent=1, default=float)
    return rep.write()


def _fmt(d):
    if abs(d) < 1e-6:
        return "$<10^{-6}$"
    m, e = ("%.1e" % d).split("e")
    return "$%s\\times10^{%d}$" % (m if d < 0 else "+" + m, int(e))


def _table(out, tab):
    L = ["reaction energies & species masses & drift per 100 steps \\\\",
         "\\hline"]
    for lab_e, lab_m, d, mx in tab:
        L.append("%s & %s & %s \\\\" % (lab_e, lab_m, _fmt(d)))
    open(os.path.join(out, "table_energy.tex"), "w").write("\n".join(L) + "\n")
    T = ["%-40s %-14s %12s %12s" % ("reaction energies", "masses",
                                     "drift/100", "max |dE/E0|")]
    for lab_e, lab_m, d, mx in tab:
        T.append("%-40s %-14s %+12.2e %12.2e" % (lab_e, lab_m, d, mx))
    open(os.path.join(out, "table_energy.txt"), "w").write("\n".join(T) + "\n")


def _fig(out, curves):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import paper_style as ps
    ps.apply()
    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    cols = {"unclosed_unconserved": ps.VERM, "unclosed_conserved": ps.ORANGE,
            "closed_unconserved": ps.PURPLE, "closed_norecomb": ps.GRAY,
            "closed_conserved": ps.BLUE, "full11": ps.BLUE}
    for name, lab_e, lab_m in ROWS:
        for i, (st, de) in enumerate(curves.get(name, [])):
            ax.plot(st, de, color=cols[name], lw=1.0,
                    ls=":" if name == "full11" else "-",
                    label="%s / %s" % (lab_e.split(" (")[0], lab_m)
                    if i == 0 else None)
    ax.axvspan(100, 300, color=ps.GRAY, alpha=0.08, lw=0)
    ax.set_xlabel("step ($2\\times10^{-11}$ s)")
    ax.set_ylabel("$E/E_0 - 1$")
    ax.legend(fontsize=5.5, loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig_energy.pdf"))
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(analyze(sys.argv[1] if len(sys.argv) > 1 else "quick"))
