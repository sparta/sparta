"""Step 8 analysis: equilibrium constants of the data files against the
NASA Glenn (CEA) database (Sec. VII F).  No SPARTA runs.

K_eq of the data files: tce_ref.Data.keq, the partition functions of
data/air.* with every mode (translation, rigid rotor with symmetry number,
harmonic vibration, electronic levels) and the reaction energy C5 of
data/air.tce.  K_eq of CEA: exp(-dG/RT) from the NASA-9 polynomials of
cea_air.thermo (standard state 1 bar), converted to number density.
Also: a Morse vibrational ladder truncated at the dissociation energy
(constants of Huber and Herzberg), and the legacy convention of using the
activation energy as the reaction energy.

Produces, in results/step8_thermo/:
  fig_keq_cea.pdf, keq_cea.txt, numbers.json, results.txt
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import pkg                 # noqa: E402
import tce_ref as R        # noqa: E402

STEP = "step8_thermo"
NA = 6.02214076e23
P0 = 1.0e5                 # Pa, standard state of the database
RXN = (("N2 + N2 --> N + N + N2", "N$_2$ $\\rightleftharpoons$ 2N"),
       ("O2 + O2 --> O + O + O2", "O$_2$ $\\rightleftharpoons$ 2O"),
       ("NO + NO --> N + O + NO", "NO $\\rightleftharpoons$ N + O"),
       ("N2 + O --> NO + N", "N$_2$ + O $\\rightleftharpoons$ NO + N"),
       ("NO + O --> O2 + N", "NO + O $\\rightleftharpoons$ O$_2$ + N"))
# Morse: E_v/k = theta v - theta x v (v + 1), theta the vibrational
# temperature of air.species, theta x = omega_e x_e hc/k (Huber and
# Herzberg 1979), levels kept while E_v increases and stays below D0
WEXE = {"N2": 14.32, "O2": 11.98, "NO": 14.075}       # cm^-1
C2 = 1.4387769                                        # hc/k, cm K
TS = (2000, 3000, 4000, 6000, 8000, 10000, 12000, 15000, 20000)


class Nasa9:
    """reader of NASA-9 records (McBride, Zehe and Gordon 2002 format)"""

    def __init__(self, path):
        L = [l for l in open(path).read().split("\n")
             if not l.startswith("!")]
        self.db = {}
        i = 0
        while i < len(L):
            l = L[i]
            if l and l[0] != " " and len(l) > 18:
                name = l[:18].split()[0]
                hdr = L[i + 1]
                n = int(hdr[:2])
                ivs = []
                for k in range(n):
                    rng = L[i + 2 + 3 * k]
                    f = lambda s: [float(s[j * 16:(j + 1) * 16].replace(
                        "D", "E")) if s[j * 16:(j + 1) * 16].strip() else 0.0
                        for j in range(5)]
                    c1, c2 = f(L[i + 3 + 3 * k]), f(L[i + 4 + 3 * k])
                    ivs.append((float(rng[:11]), float(rng[11:22]),
                                c1 + c2[:2], c2[3], c2[4], float(rng[65:80])))
                self.db[name] = dict(hf298=float(hdr[65:80]), ivs=ivs)
                i += 2 + 3 * n
            else:
                i += 1

    def g_rt(self, s, T):
        for tlo, thi, a, b1, b2, _ in self.db[s]["ivs"]:
            if tlo - 1e-6 <= T <= thi + 1e-6:
                a1, a2, a3, a4, a5, a6, a7 = a
                h = (-a1 / T ** 2 + a2 * math.log(T) / T + a3 + a4 * T / 2 +
                     a5 * T ** 2 / 3 + a6 * T ** 3 / 4 + a7 * T ** 4 / 5 +
                     b1 / T)
                s_ = (-a1 / (2 * T ** 2) - a2 / T + a3 * math.log(T) +
                      a4 * T + a5 * T ** 2 / 2 + a6 * T ** 3 / 3 +
                      a7 * T ** 4 / 4 + b2)
                return h - s_
        raise ValueError("T = %g outside the database range for %s" % (T, s))

    def kc(self, reac, prod, T):
        dg = (sum(self.g_rt(s, T) for s in prod) -
              sum(self.g_rt(s, T) for s in reac))
        dn = len(prod) - len(reac)
        return math.exp(-dg) * (P0 / (R.KB * T)) ** dn

    def dh0(self, reac, prod):
        """0-K reaction energy per molecule (J): Hf(298) - [H(298)-H(0)]"""
        h0 = lambda s: self.db[s]["hf298"] - self.db[s]["ivs"][0][5]
        return (sum(h0(s) for s in prod) - sum(h0(s) for s in reac)) / NA


def species_of(r):
    """reactants and products of the equilibrium (third body removed)"""
    if r["type"] == "D":
        return [r["reac"][0]], list(r["prod"][:2])
    return list(r["reac"]), list(r["prod"])


def q_morse(theta, thx, d0, T):
    q, v, prev = 0.0, 0, -1.0
    while True:
        e = theta * v - thx * v * (v + 1.0)
        if e <= prev or e > d0:
            return q
        q += math.exp(-e / T)
        prev, v = e, v + 1


def morse_factor(D, reac, prod, T):
    """K_eq with Morse ladders over K_eq with harmonic ladders"""
    f = 1.0
    d0 = {"N2": 1.56362e-18, "O2": 8.19608e-19, "NO": 1.04089e-18}
    for s, sgn in [(x, +1) for x in reac] + [(x, -1) for x in prod]:
        if s in WEXE:
            th = D.sp[s]["vibtemp"]
            qh = 1.0 / (1.0 - math.exp(-th / T))
            qm = q_morse(th, WEXE[s] * C2, d0[s] / R.KB, T)
            f *= (qh / qm) ** sgn
    return f


def analyze(size="quick"):
    rep = pkg.Report(STEP)
    out = os.path.join(pkg.results_dir(), STEP)
    os.makedirs(out, exist_ok=True)
    D = R.Data()
    N9 = Nasa9(os.path.join(HERE, "cea_air.thermo"))
    num = {}
    for lab, _ in RXN:
        r = D.reaction(lab)
        reac, prod = species_of(r)
        d0 = N9.dh0(reac, prod)
        dev = -r["C5"] / d0 - 1
        rep.check(abs(dev) < 5e-4, "%-24s reaction energy of data/air.tce "
                  "%.5e J, CEA 0-K value %.5e J (%+.3f%%)"
                  % (lab, -r["C5"], d0, 100 * dev),
                  meas=-r["C5"], ref=d0, sd=0.0, tol=5e-4 * abs(d0))
    Tf = np.arange(2000.0, 20000.1, 50.0)
    curves = {}
    lines = ["%-24s" % "K_eq(files)/K_eq(CEA)" + "".join("%8d" % T
                                                          for T in TS)]
    for lab, tex in RXN:
        r = D.reaction(lab)
        reac, prod = species_of(r)
        ratio = lambda T: D.keq(r, T) / N9.kc(reac, prod, T)
        v = np.array([ratio(T) for T in Tf])
        vm = v * np.array([morse_factor(D, reac, prod, T) for T in Tf])
        curves[lab] = (v, vm)
        lines.append("%-24s" % lab + "".join("%8.3f" % ratio(T) for T in TS))
        m10 = Tf <= 10000
        num[lab] = dict(at=dict((str(T), ratio(T)) for T in TS),
                        range_2_10kK=[v[m10].min(), v[m10].max()],
                        morse_factor=dict((str(T), morse_factor(
                            D, reac, prod, T)) for T in (10000, 15000, 20000)))
        rep.info("%-24s K_eq(files)/K_eq(CEA): %.3f-%.3f over 2-10 kK, %.3f "
                 "at 15 kK, %.3f at 20 kK; Morse factor %.3f/%.3f/%.3f at "
                 "10/15/20 kK" % ((lab, v[m10].min(), v[m10].max(),
                                   ratio(15000.0), ratio(20000.0)) +
                                  tuple(morse_factor(D, reac, prod, T)
                                        for T in (10000, 15000, 20000))))
    # legacy convention: reaction energy = -(activation energy)
    leg = {}
    for lab in ("N2 + O --> NO + N", "NO + O --> O2 + N"):
        r = dict(D.reaction(lab))
        r["C5"] = -r["Ea"] if lab != "N2 + O --> NO + N" else -5.175e-19
        reac, prod = species_of(r)
        leg[lab] = np.array([D.keq(r, T) / N9.kc(reac, prod, T) for T in Tf])
        b = D.keq(r, 2000.0) / D.keq(D.reaction(lab), 2000.0)
        num["legacy_" + lab] = b
        rep.info("%-24s legacy reaction energy %.4e J: K_eq x%.2f at 2000 K"
                 % (lab, -r["C5"], b))
    v = curves["N2 + N2 --> N + N + N2"]
    i20 = np.argmin(abs(Tf - 20000.0))
    share = math.log(v[0][i20] / v[1][i20]) / math.log(v[0][i20])
    num["N2_morse_share_log_20kK"] = share
    rep.info("N2 at 20 kK: the Morse ladder accounts for %.2f of the "
             "logarithmic discrepancy" % share)
    open(os.path.join(out, "keq_cea.txt"), "w").write("\n".join(lines) + "\n")
    _fig(out, Tf, curves, leg)
    with open(os.path.join(out, "numbers.json"), "w") as fh:
        json.dump(num, fh, indent=1, default=float)
    return rep.write()


def _fig(out, Tf, curves, leg):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import paper_style as ps
    ps.apply()
    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    ax.axhspan(0.95, 1.05, color=ps.GRAY, alpha=0.12, lw=0)
    cols = (ps.BLUE, ps.VERM, ps.GREEN, ps.ORANGE, ps.PURPLE)
    for (lab, tex), c in zip(RXN, cols):
        v, vm = curves[lab]
        ax.plot(Tf / 1e3, v, color=c, lw=1.2, label=tex)
        if lab.split()[0] in WEXE and "-->" in lab and lab.count("+") == 3:
            ax.plot(Tf / 1e3, vm, color=c, lw=1.0, ls=":")
        if lab in leg:
            ax.plot(Tf / 1e3, leg[lab], color=c, lw=1.0, ls="--")
    ax.set_yscale("log")
    ax.set_xlabel("$T$ (kK)")
    ax.set_ylabel("$K_{\\mathrm{eq}}^{\\mathrm{files}}/"
                  "K_{\\mathrm{eq}}^{\\mathrm{CEA}}$")
    ax.legend(fontsize=6.5, loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig_keq_cea.pdf"))
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(analyze())
