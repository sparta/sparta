"""Step 3 run list: detailed balance of the derived reverse reactions
(Secs. V and VII C).

  zeld   Zeldovich pair N2 + O <-> NO + N at 8, 10, 12, 15 kK (Fig. 7)
  rec    recombination in a dense reservoir (Table VII): atomic third
         body, K_eq fit applied to the pair (keq_file), molecular third
         body, molecular third body with continuous vibration, at REC_T;
         and the atomic case at 6765 K, where recombination dominates
  mech   11-species reservoir with the full mechanism at 15 and 8 kK
         (the 28 dissociation and exchange pairs; new figure)
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pkg import Run, DATA, size_params  # noqa: E402

STEP = "step3_balance"
# NSTEP: Zeldovich reservoirs (dilute); NSTR: recombination reservoirs
# (dense, about 5 collisions per particle per step); NSTM: 11-species
# reservoirs (dense)
# REC_T: temperature of the recombination configurations.  The paper uses
# 15 kK, where recombination tallies are about 1/2400 of the dissociation
# tallies at this density, too few for a short run; at 10 kK they are
# 1/40.  Lower temperatures resolve the ratio sooner but make it more
# sensitive to the temperature (d ln K_eq/dT = D0/kT^2: 0.2% per 10 K at
# 10 kK), so every reference is evaluated at the measured mean
# temperature of the run.
SIZES = {"quick": dict(seeds=2, NP=40000, NSTEP=4000, NSTR=6000,
                       REC_T=10000.0, NPM=40000, NSTM=400)}
FILES = ("zeldovich.tce", "recomb_atom.tce", "recomb_atom_fwd.tce",
         "recomb_mol.tce", "park_recomb.keq")
NLOW, NHIGH = 7.07043e22, 1.414086e25
ZELD_T = (8000.0, 10000.0, 12000.0, 15000.0)
REC = {  # name: (reaction file, VIB, react_modify extra, T; None = REC_T)
    "atom": ("recomb_atom.tce", "discrete", "recomb yes", None),
    "keq": ("recomb_atom_fwd.tce", "discrete",
            "reverse auto keq_file park_recomb.keq", None),
    "mol": ("recomb_mol.tce", "discrete", "recomb yes", None),
    "molsmooth": ("recomb_mol.tce", "smooth", "recomb yes", None),
    "atom6765": ("recomb_atom.tce", "discrete", "recomb yes", 6765.0),
}


def runs(size="quick"):
    sz = size_params(STEP, SIZES, size)
    out = []
    for s in range(sz["seeds"]):
        seed = 32345 + 1000 * s
        for T in ZELD_T:
            v = dict(DATA=DATA, SEED=seed, T=T, NRHO=NLOW, NP=sz["NP"],
                     NSTEP=sz["NSTEP"], TCE="zeldovich.tce", VIB="discrete",
                     RMOD="recomb yes")
            out.append(Run(STEP, "zeld_T%d_s%d" % (T, s), "in.zeld", v,
                           files=FILES))
        for name, (tce, vib, rmod, T) in REC.items():
            T = T or sz["REC_T"]
            v = dict(DATA=DATA, SEED=seed, T=T, NRHO=NHIGH, NP=sz["NP"],
                     NSTEP=sz["NSTR"], TCE=tce, VIB=vib, RMOD=rmod)
            out.append(Run(STEP, "rec_%s_s%d" % (name, s), "in.res", v,
                           files=FILES))
        for T in (15000.0, 8000.0):
            v = dict(DATA=DATA, SEED=seed, T=T, NRHO=NHIGH, NP=sz["NPM"],
                     NSTEP=sz["NSTM"])
            out.append(Run(STEP, "mech_T%d_s%d" % (T, s), "in.mech", v))
    return out
