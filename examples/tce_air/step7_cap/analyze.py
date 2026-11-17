"""Step 7 analysis: the probability cap with reverse reactions.

Produces, in results/step7_cap/:
  table10_captest.tex/.txt  Table X: realized equilibrium constant of
                            NO <-> N + O relative to the partition-function
                            value, per reactor configuration
  numbers.json, results.txt

The frozen-reservoir rates vary from seed to seed by several times the
counting error of their tallies (one 20 000-particle reservoir at 8000 K
scatters by about 5% against a 0.7% counting error, at 10-15 kK by
1-2%), because each reservoir keeps features of its initial random
sample that carry a large share of these rates:
  - spin classes.  Non-reacting collisions conserve spin (air.elec sets
    no exception), and a frozen reservoir performs no reaction, so the
    number of NO quartets (a4Pi, b4Sigma-), O singlets and N doublets
    stays at its initial sample for the whole run.  At 8000 K the 0.2%
    of NO in quartet states carries 13% of NO + O and 10% of NO + N
    dissociation, the 3% of O in singlets 18% of NO + O, the 8% of N in
    doublets 50% of NO + N (tce_ref, capped average by reactant state);
    a reservoir holds about 14 quartet NO molecules.
  - vibration.  NO relaxes its vibration once in about 1400 collisions,
    a few times per run, so the hot tail of its initial sample persists.
The check takes the error of the mean from the spread over seeds, and
allows 2% beyond it: the measurements run 0.7-2.7% below the references
at 8-15 kK for both reactions.  The relaxation model does not bias the
rates: in 100 000-particle reservoirs at 8000 K every mode stays at the
Boltzmann distribution of the measured temperature, and the rates agree
with the capped fraction to 1-2% whether NO relaxes at the rates of the
data or at every collision (rotation, vibration or both).

References (tce_ref.py, from the data files only):
  frozen reservoirs: the capped fraction <min(P,1)>/<P> of the forward
    microcanonical probability (ExchangeReverse with forward=True, which
    is calibrated to the Arrhenius fit), at the measured temperature; the
    recombination is built from the capped forward numerator (Eq. capnum),
    so its rate relative to k_fit/K_eq is the same fraction
  reactors: the temperature at which the mean composition and the
    conserved total energy (thermal energy of all modes plus formation
    energy) agree; K_eq(T) of the partition functions
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
from runs import runs, STEP, FRZ_T  # noqa: E402

V = (2.75e-6) ** 3
FWD = ("NO + O --> N + O + O", "NO + N --> N + O + N")
REV = {"NO + O --> N + O + O": "N + O --> NO + O",
       "NO + N --> N + O + N": "N + O --> NO + N"}
ROWS = (("precap", "forward-only cap (build before the cap correction)"),
        ("cap", "cap correction, Eq. (capnum)"),
        ("lown", "cap correction, a tenth of the density"),
        ("slow", "rates divided by 30 (no cap)"),
        ("fast", "cap correction, NO relaxation at every collision"))


def check_files(rep):
    air = [" ".join(l.split()) for l in open(os.path.join(R.DATA, "air.tce"))]
    lines = [" ".join(l.split()) for l in open(os.path.join(HERE, "cap.tce"))
             if l.strip() and not l.startswith("#")]
    rep.check(all(l in air for l in lines), "cap.tce agrees with data/air.tce")
    a = [l.split() for l in open(os.path.join(HERE, "nofast.species"))
         if l.strip() and not l.startswith("#")]
    b = [l.split() for l in open(os.path.join(R.DATA, "air.species"))
         if l.strip() and not l.lstrip().startswith("#") and
         len(l.split()) >= 10]
    same = len(a) == len(b) and all(
        x[0] == y[0] and all(float(p) == float(q) for k, (p, q) in
                             enumerate(zip(x[1:], y[1:]))
                             if not (x[0] == "NO" and k + 1 in (4, 6)))
        for x, y in zip(a, b))
    rep.check(same, "nofast.species equals data/air.species except the "
              "relaxation probabilities of NO")


def stats_array(run):
    head, rows = pkg.stats(run.log)
    a = np.array(rows)
    return {h: a[:, i] for i, h in enumerate(head)}


def equilibrium_T(D, N, E):
    """T at which sum_s N_s (e_mean(s, T) + E_f,s) = E"""
    f = lambda T: sum(N[s] * (D.e_mean(s, T) + R.EFORM[s]) for s in N) - E
    lo, hi = 1000.0, 50000.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def analyze(size="quick"):
    rep = pkg.Report(STEP + " (" + size + ")")
    out = os.path.join(pkg.results_dir(), STEP)
    os.makedirs(out, exist_ok=True)
    check_files(rep)
    D = R.Data()
    rl = runs(size)
    num = {}

    # ---- frozen reservoirs
    refs = {lab: R.ExchangeReverse(D, D.reaction(lab), forward=True)
            for lab in FWD}
    for T in FRZ_T:
        sel = [q for q in rl if q.tag.startswith("frz_T%d_s" % T)
               and q.done()]
        if not sel:
            continue
        # per seed and reaction: realized rate over the reference, its
        # counting error, the measured fraction and the capped fraction
        fwd = dict((lab, []) for lab in FWD)
        rec = dict((lab, []) for lab in FWD)
        Tms = []
        for r in sel:
            st = stats_array(r)
            Ts = st["c_temp"][1:]               # every 20 steps
            Tm = float(np.mean(Ts))
            Tms.append(Tm)
            t = r.variables["NSTEP"] * r.variables["DT"]
            fnum = r.variables["NRHO"] * V / r.variables["NP"]
            n = r.variables["NRHO"]
            tal = pkg.tallies(r.log)
            # composition actually created (mole fractions of the sample)
            x = {"O": st["c_spc[1]"][0], "N": st["c_spc[2]"][0],
                 "NO": st["c_spc[3]"][0]}
            tot = sum(x.values())
            x = {k: v / tot for k, v in x.items()}
            for lab in FWD:
                rx = D.reaction(lab)
                M = rx["reac"][1]
                capf = refs[lab].mean_p(Tm) / refs[lab].mean_p(Tm, cap=False)
                # the tallies integrate the rate over the run, and the
                # rates vary as T^5-T^10: average the fit over the
                # temperature series (capf varies slowly, at Tm)
                F = tal.get(lab, 0)
                if T >= 8000.0 and F > 0:
                    k = F * fnum / (x["NO"] * n * x[M] * n * V * t)
                    kfit = float(np.mean([R.arrhenius(rx, q) for q in Ts]))
                    fwd[lab].append((k / kfit / capf, 1 / math.sqrt(F),
                                     k / kfit, capf))
                B = tal.get(REV[lab], 0)
                if T < 8000.0 and B > 0:
                    k = B * fnum / (x["N"] * n * x["O"] * n * x[M] * n * V * t)
                    kdb = float(np.mean([R.arrhenius(rx, q) / D.keq(rx, q)
                                         for q in Ts]))
                    rec[lab].append((k / kdb / capf, 1 / math.sqrt(B),
                                     k / kdb, capf))
        Tm = float(np.mean(Tms))
        for lab in FWD:
            for kind, vals in (("fwd", fwd[lab]), ("rec", rec[lab])):
                if not vals:
                    continue
                y = [v[0] for v in vals]
                m = float(np.mean(y))
                # error of the mean: counting statistics, or the spread
                # over seeds if larger (it is, by several times)
                e = math.sqrt(sum((v[0] * v[1]) ** 2 for v in vals)) / len(y)
                if len(y) > 1:
                    e = max(e, float(np.std(y, ddof=1)) / math.sqrt(len(y)))
                frac = float(np.mean([v[2] for v in vals]))
                capf = float(np.mean([v[3] for v in vals]))
                num["frz_%s_%s_%d" % (kind, lab.split()[2], T)] = dict(
                    ratio=[m, 2 * e], fraction=frac, capped_fraction=capf,
                    T=Tm, seeds=y)
                what = ("%s at %5.0f K: %.3f of the fit" % (lab, Tm, frac)
                        if kind == "fwd" else
                        "%s at %5.0f K: %.3f of k_fit/K_eq" % (REV[lab], Tm,
                                                              frac))
                rep.check(abs(m - 1) <= 3 * e + 0.02,
                          "%s, capped %sfraction %.3f: ratio %.3f +/- %.3f "
                          "(%d seeds)" % (what, "forward " if kind == "rec"
                                          else "", capf, m, 2 * e, len(y)))

    # ---- reactors
    tab = []
    for key, lab in ROWS:
        sel = [r for r in rl if r.tag.startswith("rx_%s_s" % key)
               and r.done()]
        if not sel:
            if key == "precap":
                rep.info("forward-only cap: not run (needs the 'precap' "
                         "build, see README)")
            continue
        ratios, errs = [], []
        for r in sel:
            st = stats_array(r)
            steps = st["Step"]
            keep = np.where(steps >= 0.5 * steps[-1])[0]
            fnum = r.variables["NRHO"] * V / r.variables["NP"]
            E = st["v_etot"][0]
            drift = abs(st["v_etot"][-1] / E - 1)
            # the ratio at every snapshot of the second half, from the
            # composition and the temperature that the conserved energy
            # implies for it: near equilibrium the composition fluctuates,
            # and one net dissociation in 10 000 particles moves the ratio
            # by about 1% (through the temperature it implies); the error
            # of the mean follows from the autocorrelation time of the
            # series, which is long when the reactions are slow
            rs, Ts = [], []
            for i in keep:
                N = {"O": st["c_spc[1]"][i], "N": st["c_spc[2]"][i],
                     "NO": st["c_spc[3]"][i]}
                T = equilibrium_T(D, N, E)
                K = (N["N"] * N["O"] / N["NO"]) * fnum / V
                rs.append(K / D.keq(D.reaction("NO + O --> N + O + O"), T))
                Ts.append(T)
            ratio, err = pkg.mean_err_corr(rs)
            ratios.append(ratio)
            errs.append(err)
            num["rx_%s" % r.tag[3:]] = dict(
                ratio=ratio, err=err, T=float(np.mean(Ts)),
                energy_drift=drift,
                T_tr=float(st["c_temp"][keep].mean()))
            rep.info("%-18s T %.0f K (translational %.0f K), K_eq realized/"
                     "partition function %.3f +/- %.3f, energy drift %.1e"
                     % (r.tag, np.mean(Ts), st["c_temp"][keep].mean(), ratio,
                        2 * err, drift))
        m = float(np.mean(ratios))
        e = math.sqrt(sum(x * x for x in errs)) / len(errs)
        if len(ratios) > 1:
            e = max(e, float(np.std(ratios, ddof=1)) / math.sqrt(len(ratios)))
        tab.append((lab, m, ratios))
        if key in ("cap", "lown", "slow", "fast"):
            rep.check(abs(m - 1) <= 3 * e + 0.02,
                      "%s: realized/partition-function K_eq %.3f +/- %.3f "
                      "(2 sd)" % (lab, m, 2 * e))
    _table(out, tab)
    with open(os.path.join(out, "numbers.json"), "w") as fh:
        json.dump(num, fh, indent=1, default=float)
    return rep.write()


def _table(out, tab):
    L = ["configuration & realized/partition-function \\\\", "\\hline"]
    T = ["%-55s %s" % ("configuration", "realized/partition function")]
    for lab, m, rs in tab:
        s = "%.2f" % m if len(rs) == 1 else "%.2f (%s)" % (
            m, ", ".join("%.2f" % x for x in rs))
        L.append("%s & %s \\\\" % (lab.replace("(capnum)", "(\\ref{eq:capnum})"), s))
        T.append("%-55s %s" % (lab, s))
    open(os.path.join(out, "table10_captest.tex"), "w").write("\n".join(L) + "\n")
    open(os.path.join(out, "table10_captest.txt"), "w").write("\n".join(T) + "\n")


if __name__ == "__main__":
    sys.exit(analyze(sys.argv[1] if len(sys.argv) > 1 else "quick"))
