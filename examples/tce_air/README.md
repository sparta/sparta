# tce_air: the paper's claims, step by step

This directory checks the main claims of the paper *A rate-preserving
microcanonical total collision energy model for discrete internal states
and detailed-balance reverse reactions* with SPARTA runs on an ordinary
computer.  Each step runs a few small simulations of air chemistry and
compares what SPARTA does with a reference that is computed independently,
from the data files and the formulas of the paper, by `tce_ref.py`.  Every
comparison is printed as a PASS or FAIL line.

It is meant for reviewing the code and the physics.  You need no knowledge
of SPARTA internals: each step is one directory with its input decks
(`in.*`), its run list (`runs.py`) and its analysis (`analyze.py`), and
the analysis says in words what it compares.

## What you need

- A Linux or macOS machine with 4 cores (more is faster), an MPI C++
  compiler (`mpicxx`), `make`, and Python 3 with `numpy`, `scipy` and
  `matplotlib`.
- About 2 hours for all steps with 4 cores (about 7 core-hours; the ten
  recombination runs of `step3_balance`, half an hour each, are three
  quarters of it).

## Quick start

```bash
cd examples/tce_air
./build.sh                       # builds bin/spa_release from this branch
./tce_air.py run all --jobs 4 --np 1
./tce_air.py analyze all
```

`run` skips runs that already finished, so it can be repeated after an
interruption.  `analyze` writes everything to `results/<step>/` and
returns the number of failed checks.  On a machine with more cores, run
more jobs at a time: `--jobs 6 --np 2` on 12 cores.

Two optional builds add comparisons with earlier treatments:

```bash
./build.sh --compare             # bin/spa_compare: the earlier TCE treatments
./build.sh --precap              # bin/spa_precap: before the cap correction
```

Without them the runs that need them are skipped, and their checks are
left out.

## What each step checks

Each step corresponds to a part of the paper.  The step's `runs.py` lists
its runs, `analyze.py` its checks; `results/<step>/results.txt` has the
outcome.

| step | paper | what is checked |
|---|---|---|
| `step1_rates` | Secs. III-IV, VII A; Figs. 1, 3; Table VI | The equilibrium dissociation rate of O2 + N2 reproduces the input Arrhenius rate with the microcanonical factor, for vibrational, electronic and joint ladders and discrete rotation; the earlier treatments (adding electronic energy to the collision energy, leaving it out, instantaneous vibrational DOF) miss it by the amounts the paper derives. |
| `step2_noneq` | Sec. VII B; Figs. 4-6; Table VII | With vibrational or electronic populations frozen away from the translational temperature, the rate follows an independent quadrature of the ladder average; comparison with Bird's QK model, Park's two-temperature model and QCT data for N2 + N2. |
| `step3_balance` | Secs. V, VII C; Fig. 7; Table VIII | Derived reverse reactions balance their forward partners: the Zeldovich exchange, three-body recombination with atomic and molecular third bodies, with an equilibrium-constant fit, and the full 11-species mechanism. |
| `step4_box` | Sec. VII C; Fig. 8 | A closed reacting box with the Zeldovich pair relaxes to the equilibrium composition from both sides, thermostatted and adiabatic, and the internal modes of every species stay at the translational temperature. |
| `step5_energy` | Sec. VII D; Table IX | Total energy is conserved through dense reaction cycling only when reaction energies close around cycles and species masses are conserved; each defect alone makes the energy drift. |
| `step6_noneq_reverse` | Sec. VII E; Table X | The reverse rate responds to the internal temperature with no cell temperature in the model, mode by mode, and recombination responds to the electronic populations, as the microcanonical tables predict. |
| `step7_cap` | Secs. IV D, VII G; Table XI | Where fitted rates exceed the collision rate, the capped forward rate and its reverse follow the capped fraction of the microcanonical probability, and a reacting gas still reaches the partition-function equilibrium. |
| `step8_thermo` | Sec. VII F; Fig. 9 | The equilibrium constants of the data files against the NASA Glenn (CEA) database; no SPARTA runs. |

Figure and table numbers are those of the draft of 2026-11-24.

## Figures and tables

Each step writes its figures and tables to `results/<step>/`, named by
content.  Every step except `step5_energy` also writes `summary.pdf`: one
row per check, the measured value over its reference with the 2-sd bar,
against the pass window (gray).  `checks.json` has the same numbers.

| step | file | in the paper |
|---|---|---|
| `step1_rates` | `fig_progression.pdf` | Fig. 1 |
| | `fig_eq_discretevib.pdf` | Fig. 3 |
| | `table_reservoir.tex/.txt` | Table VI |
| `step2_noneq` | `fig_noneq.pdf` | Fig. 4 |
| | `fig_reservoir.pdf` | Fig. 5 |
| | `fig_qk_compare.pdf` | Fig. 6 |
| | `table_qct.tex/.txt` | Table VII |
| `step3_balance` | `fig_zeldovich.pdf` | Fig. 7 |
| | `table_recomb.tex/.txt`, `fig_recomb.pdf` | Table VIII (and as a figure) |
| | `fig_mechanism.pdf`, `mechanism_*K.txt` | the 11-species reservoir (supplement) |
| `step4_box` | `fig_closedbox.pdf` | Fig. 8 |
| `step5_energy` | `table_energy.tex/.txt`, `fig_energy.pdf` | Table IX (and as a figure) |
| `step6_noneq_reverse` | `table_cellT.tex/.txt`, `fig_cellT.pdf` | Table X (and as a figure) |
| `step7_cap` | `table_captest.tex/.txt`, `fig_cap.pdf` | Table XI; the frozen-reservoir numbers of Sec. VII G |
| `step8_thermo` | `fig_keq_cea.pdf` | Fig. 9 |

The figures are made from the results already on disk:
`./tce_air.py analyze all` redraws them without running SPARTA, so a
`results/` directory from an earlier version of the package gets the new
figures this way.  Earlier versions numbered the files after an older
draft (`fig03_progression.pdf`, `table7_recomb.tex`, ...); the analysis
removes those.

## Reading the results

`results/<step>/results.txt` lists every check:

```
  [PASS] TV=  300 K microcanonical: 0.538 +/- 0.012, analytic 0.538
```

The measured value carries two standard deviations of its counting
statistics (or of the spread over seeds, if larger); a check passes when
the measurement agrees with the reference within three standard
deviations plus a small allowance stated in the analysis.  Lines without
PASS or FAIL are information (numbers the paper quotes, or values a check
does not need).  Each step also writes its figures (`fig*.pdf`), tables
(`table*.tex`, `.txt`) and `numbers.json` (see "Figures and tables").

The quick size of the package resolves the claims at the few-percent
level.  The paper's numbers come from the same decks at larger sizes,
run on GPUs (the `paper-runs` branch of the paper's repository).

## What to expect

With this branch at the quick size all checks pass.  Three things are
worth knowing:

- `step7_cap`: the forward rates of NO + O and NO + N dissociation, whose
  Arrhenius fits exceed the collision rate so that most reacting
  collisions are capped at probability one, run 0.7-2.7% below the
  capped fraction computed by `tce_ref.py` at 8000-15 000 K; the check
  allows 2%.  Single frozen reservoirs of this step scatter by several
  times their counting error (about 5% at 8000 K).  A frozen reservoir
  performs no reaction and its collisions conserve spin, so it keeps its
  initial random numbers of excited spin states (about 14 NO molecules
  in quartet states, which carry 10-13% of the rate at 8000 K), and NO
  relaxes its vibration only a few times per run.  The step therefore
  runs six seeds per temperature and takes the error from their
  spread.  The reactors of the same step reach the
  partition-function equilibrium to 1-4%.
- `step3_balance`, the 11-species reservoir: every dissociation and
  exchange pair with enough tallies agrees with detailed balance within
  3 standard deviations (plus 1%), with the composition and temperature
  that each sample actually has.  The three associative-ionization pairs
  (O + N -> NO+ + e and the like) do not: their tally ratios are 13-27%
  above the Saha prediction at 15 000 K and about 80% at 8000 K.  Their
  reverse, dissociative recombination, is capped at low electron energy,
  a loss the forward reaction does not see; the paper discusses this,
  and the analysis reports these pairs without checking them.
- The recombination checks of `step3_balance` have the largest
  statistical errors of the package at the quick size: 5-9% at two
  standard deviations (42% at 6765 K, where few dissociations occur).

Two limitations are known and documented in the paper:

- A barrierless reverse reaction (here NO + N -> N2 + O) has
  probabilities above one at low collision energy; their cap lowers the
  reverse rate slightly (0.6% at 8000 K), and this is not carried back
  into the forward reaction.  `step3_balance` checks against the capped
  expectation.
- Dissociation produces three particles, and SPARTA scatters them with
  the usual DSMC convention (two of them leave together), which is not
  the state distribution assumed by the derived recombination.  In a
  dense reactor where only dissociation and recombination act, the
  realized equilibrium constant then depends on the number of particles
  per cell (a few percent at 1000-5000 per cell); the GPU runs measure
  this.

## How the references are computed

`tce_ref.py` reads the data files (`data/air.*`) and implements the
paper's formulas directly: partition functions and equilibrium
constants, the TCE probability and its microcanonical ladder average, the
derived reverse probability with its calibration and the cap at one, and
the three-body recombination average.  It does not call SPARTA or use any
of its tables.  Where the paper says "quadrature", `tce_ref.py` computes
the integral numerically; each reference function also checks itself
against a closed form where one exists (for example, the calibration
constant of a reverse table must not depend on temperature).

## Files

```
build.sh        builds bin/spa_release (and optionally spa_compare, spa_precap)
tce_air.py      list, run and analyze the steps
pkg.py          runs, logs and reports shared by the steps
tce_ref.py      the independent reference calculations
paper_style.py  figure style
stepN_*/        one step: in.* decks, runs.py, analyze.py, data files
results/        written by run and analyze
```

Each run writes `results/<step>/<run>/log.sparta` and `run.json` (how the
run was made); a run whose definition changes is run again.
