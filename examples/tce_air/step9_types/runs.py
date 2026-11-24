"""Step 9 run list: reaction types beyond dissociation.

One frozen reservoir per forward/reverse pair of data/air.tce, with the
four species of the pair present, at T in TYPES_T:
  exch   exchange               N2 + O --> NO + N   (reverse NO + N --> N2 + O)
  ai     associative ionization O + N --> NO+ + e   (reverse: dissociative
                                                     recombination NO+ + e)
  cx     charge exchange        O+ + N2 --> O + N2+ (reverse N2+ + O)

Treatments of the forward reaction:
  micro  this branch, microcanonical factor with the discrete vibrational
         and electronic ladders (release build, in.types); the derived
         reverse is tallied in the same runs
  incl   electronic energy added to the collision energy, vibration left
         out (comparison build, in.types_cmp with VIB no)
  inst   Bird's instantaneous vibrational DOF, no electronic states
         (comparison build, in.types_noelec_cmp); only for pairs with a
         molecule

Electron-impact ionization (N + e --> N+ + e + e) is not here: data/air.tce
leaves it out because its fits have eta = -3.8, below the TCE bound
eta > -3/2 of an atom-electron pair (see the comment in data/air.tce).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pkg import Run as _Run, DATA, size_params  # noqa: E402

STEP = "step9_types"
# NSTEP: the forward tallies at 8000 K are the fewest (about 400 per run of
# 20 000 steps for charge exchange)
SIZES = {"quick": dict(seeds=2, NP=40000, NSTEP=20000)}
NRHO = 7.07043e22
TYPES_T = (8000.0, 10000.0, 15000.0, 20000.0)

# pair: (species S1-S4, fractions X1-X4, forward, reverse, label)
PAIRS = {
    "exch": (("N2", "O", "NO", "N"), (0.25, 0.25, 0.25, 0.25),
             "N2 + O --> NO + N", "NO + N --> N2 + O", "exchange"),
    # few electrons: their collisions dominate the cost of the run, and
    # the dissociative recombination still tallies by the ten thousand
    "ai": (("O", "N", "NO+", "e"), (0.45, 0.45, 0.05, 0.05),
           "O + N --> NO+ + e", "NO+ + e --> O + N",
           "associative ionization"),
    "cx": (("O+", "N2", "O", "N2+"), (0.25, 0.25, 0.25, 0.25),
           "O+ + N2 --> O + N2+", "N2+ + O --> O+ + N2", "charge exchange"),
}

TREAT = {
    # name: (executable, deck, react_modify of the comparison build, VIB)
    "micro": ("release", "in.types", None, None),
    "incl": ("compare", "in.types_cmp",
             "partial_energy no vib_energy dof elec_energy yes", "no"),
    "inst": ("compare", "in.types_noelec_cmp",
             "partial_energy no vib_energy dof", None),
}

FILES = ("exch.tce", "exch_fwd.tce", "ai.tce", "ai_fwd.tce", "cx.tce",
         "cx_fwd.tce")


def has_molecule(pair):
    return any(s in ("N2", "O2", "NO", "N2+", "O2+", "NO+")
               for s in PAIRS[pair][0][:2])


def treatments(pair):
    return [t for t in TREAT if t != "inst" or has_molecule(pair)]


def runs(size="quick"):
    sz = size_params(STEP, SIZES, size)
    out = []
    for s in range(sz["seeds"]):
        seed = 52345 + 1000 * s
        for pair, (sp, x, fwd, rev, lab) in PAIRS.items():
            for tr in treatments(pair):
                kind, deck, rmod, vib = TREAT[tr]
                for T in TYPES_T:
                    v = dict(DATA=DATA, SEED=seed, T=T, NRHO=NRHO,
                             NP=sz["NP"], NSTEP=sz["NSTEP"],
                             TCE=pair + (".tce" if tr == "micro"
                                         else "_fwd.tce"))
                    for i in range(4):
                        v["S%d" % (i + 1)] = sp[i]
                        v["X%d" % (i + 1)] = x[i]
                    if rmod:
                        v["RMOD"] = rmod
                    if vib:
                        v["VIB"] = vib
                    out.append(_Run(STEP, "%s_%s_T%d_s%d" % (pair, tr, T, s),
                                    deck, v, kind=kind, files=FILES))
    return out
