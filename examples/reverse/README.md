# Detailed-balance reverse reactions (issue #472)

This example demonstrates **backward reaction rates derived from forward
rates by detailed balance**, instead of requiring an independently-fitted
Arrhenius expression for each reverse reaction.  It implements the
enhancement requested in issue #472, following Bird (1994) Sec. 6.6 and
Boyd & Schwartzentruber (2017) Secs. 7.5.2-7.5.3.

## Why

In the TCE model every reaction - forward or reverse - is normally an
independent entry in the reaction file with its own Arrhenius triple
`(A, b, Ea)`.  A single Arrhenius form cannot exactly reproduce
`k_b(T) = k_f(T) / K_eq(T)`, because `K_eq` carries the non-Arrhenius
temperature dependence of the partition-function ratios; independently
fitted forward/backward rates therefore do not satisfy detailed balance,
and a reacting gas does not relax to the correct chemical equilibrium
composition.

## The B reaction style

A reaction with style **`B`** in the reaction file gets its rate derived
at run time from its forward partner (see `doc/src/react.rst` for the full
description).  Two kinds of pairs are supported:

- **exchange <-> exchange**: the forward partner is the exchange reaction
  with reactants and products swapped
- **recombination <-> dissociation**: `A + B -> AB + M` pairs with the
  per-partner dissociation `AB + M -> A + B + M` with the same explicit
  third body `M`

At initialization the backward activation energy (`Ea_B = Ea_F + dHf`),
reaction energy (`-dHf`), and effective DOF are seeded from the forward
line; the five numeric coefficients on a `B` line are placeholders.

An **exchange** pair is then given a temperature-free detailed-balance
table: the ratio of the two channels' microcanonical densities of
states at matched total collision energy (built with the same
machinery as the microcanonical TCE energy factors), so each
collision's backward probability follows from the forward reaction by
energy-resolved microscopic reversibility.  Its thermal average
reproduces `k_b(T) = k_f(T)/K_eq(T)` at every temperature
simultaneously, and no temperature is evaluated at run time.

A **recombination** pair gets the same treatment extended to three
bodies: its probability is resolved in the total available energy

```
w = u_pair + eps_3 + e_int,3
```

(the collision energy of the recombining pair, plus the third
particle's translational energy relative to the pair's center of mass
- which is exactly the relative translational energy of the forward
dissociation collision - plus the third particle's internal energies).
The probability is the calibrated ratio of the forward channel's
microcanonical reaction numerator at w to the density of states of
the (pair, third-body) energy decomposition, so an energetic third
body enhances recombination exactly as microscopic reversibility of
the forward dissociation demands.  As for exchange, the thermal
average reproduces `k_b(T) = k_f(T)/K_eq(T)` at every temperature and
no cell temperature is used anywhere.

Partition functions (used to calibrate the tables at initialization)
include translational, rigid-rotor rotational (with the symmetry
number sigma from the species `rotfile`), harmonic-oscillator
vibrational, and electronic (from the species `elecfile` ladder)
factors.

Reverse reactions require `react tce` with the microcanonical
`react_modify partial_energy no` coupling (otherwise B lines are
switched off with a warning), a `rotfile` entry for every rotating
species they involve, and are not available for
ionization, whose reverse rate depends on the electron temperature and
must be supplied explicitly.  The KOKKOS `tce/kk` style produces
bit-for-bit identical results (SPARTA_KOKKOS_EXACT).

## Files

- `in.reverse` - demonstration deck: thermal air box (N2/O/NO/N) at
  20000 K with all internal modes discrete, running the two pairs in
  `rev.tce`
- `rev.tce` - exchange pair `N2 + O <-> NO + N` and
  dissociation/recombination pair `N2 + N <-> N + N + N`; forward
  parameters from `data/air.tce`
- `rev_exch.tce` - exchange pair only (used by the relaxation check)
- `in.reverse_rate` - frozen-composition reservoir
  (`compute_chem_rates yes`) for measuring rates; variables `T`, `RB`,
  `NRHO`, `FNUM`
- `in.reverse_eq` - closed reacting box for the equilibrium-relaxation
  check; variables `T`, `FN2`, `FO`, `FNO`, `FN`
- `air.rot` - rotational temperatures + symmetry numbers (N2: 2.88 K,
  sigma 2; O2: 2.07 K, sigma 2; NO: 2.44 K, sigma 1)
- `air.elec` - NIST-based low-lying electronic levels for N2, O, N, NO
- `validate_reverse.py` - quantitative validation battery (below)

## Running

```
../../src/spa_serial -in in.reverse
```

Both exchange directions fire in the tallies; the exothermic reverse is
more frequent than the endothermic forward, as expected:

```
reaction N2 + O --> NO + N: 197
reaction NO + N --> N2 + O: 506
reaction N2 + N --> N + N + N: 1004
```

Three-body recombination events are too rare to appear at this density
in a 1000-step demonstration (the rate scales as density cubed); the
validation battery measures them in a 200x denser reservoir.

## Validation

`validate_reverse.py` recomputes the partition functions and equilibrium
constants independently of the SPARTA implementation (from the same data
files) and checks the measured rates against them:

```
python3 validate_reverse.py --exe ../../src/spa_serial
```

Add `--exe2 <kokkos-binary> --exe2-args "-k on -sf kk -pk kokkos react/retry
yes"` to also run the CPU/KOKKOS parity check (check 12).

1. **Exchange detailed balance**: the forward/backward tally ratio from
   frozen reservoirs matches the analytic `K_eq(T)` to 3.5% at 15000 K,
   3.6% at 10000 K, and 9.5% at 8000 K (within the tally statistics).
2. **Literature comparison**: the derived backward rate for
   `NO + N -> N2 + O` lies within a factor of 1.33-1.39 of the
   independently fitted literature rate for the same reaction in
   `data/air.tce` over 8000-15000 K - i.e. the derived rate agrees with
   the fitted one to well within the scatter of published rate models.
   The residual is the fit pair's own inconsistency: the equilibrium
   constant implied by the two published fits is x1.38-1.40 larger than
   the statistical-mechanics K_eq, so a thermodynamically consistent
   backward rate cannot (and should not) match the fitted one exactly.
3. **Recombination detailed balance**: in a dense reservoir the
   dissociation/recombination tally ratio times the atom number density
   matches the analytic volumetric `K_eq` to 1.4% (2.9% statistics) -
   with the fully microcanonical 3-body probability, no cell
   temperature involved.
4. **Equilibrium relaxation**: a closed reacting box initialized on the
   pure-reactant side and on the pure-product side relaxes toward the
   same analytic equilibrium composition (NO fraction 0.411) from both
   directions.
5. **Input sanity**: the electronic ground-state degeneracies and
   symmetry numbers that dominate the `K_eq` prefactor are present in the
   data files (for this reaction they contribute factors of ~1.8 and 2).
6. **Error paths**: a `B` reaction without a forward partner, a `B`
   recombination with a wildcard third body, and a `B` reaction under a
   QK style all abort with the intended error messages.
7. **Auto-generation**: `react_modify reverse auto` on a forward-only
   reaction file generates the same two reverse reactions as the
   hand-written `B` lines and holds the same detailed balance.
8. **External K_eq**: `react_modify keq_file` with a Park-form fit of
   the equilibrium constant implied by the published forward/reverse
   pair reproduces the published reverse rate to x1.02 - the option to
   use when a chemistry set's backward rates must be matched exactly
   rather than derived from statistical mechanics.  The matched reverse
   keeps its energy-resolved detailed-balance table and applies only the
   residual thermal factor `R(T) = K_eq_statmech/K_eq_fit`, so the
   collision-energy selectivity stays microscopically reversible and only
   the (small) K_eq discrepancy is thermal.
9. **External K_eq for recombination**: a Park fit of the analytic
   volumetric dissociation `K_eq` (1/m^3), fed via `keq_file`, reproduces
   that `K_eq` as the three-body recombination backward rate to within
   statistics - the m^6/s recombination analogue of check 8, exercising
   the volumetric-unit path the dimensionless exchange check cannot.
10. **No spurious bounds warning**: the ubiquitous `eta = -3/2`
    dissociation with one rotor sits exactly on the low-energy trend
    bound of the TCE probability, where the threshold factor is finite;
    `check_tce_bounds` must not warn that its rate is erroneous.
11. **Detailed balance under `vibrate smooth`**: with classical
    (continuum) vibration the detailed-balance table and its calibration
    target share the same vibrational temperature dependence, so the
    table does not report drift and the exchange reverse reproduces the
    classical-vibration `K_eq`.
12. **KOKKOS parity of the external-K_eq path** (only with `--exe2`):
    the cell-temperature-dependent `keq_file` residual factor produces
    reaction tallies bit-for-bit identical between the CPU and KOKKOS
    styles.
13. **Molecular third body**: the recombination `N + N -> N2 + N2` (data
    in `rev_mol.tce`) has a molecular third body N2, whose discrete
    vibrational and electronic ladders fold into the 3-body density of
    states and whose continuum rotation is a flat measure variable - the
    general case of the microcanonical recombination that check 3, with an
    atomic third body N, does not exercise.  The dissociation/recombination
    tally ratio matches the volumetric `K_eq` to within statistics.
14. **Non-equilibrium reverse rate**: the whole point of the temperature-free
    per-collision construction is far-from-equilibrium behaviour.  The
    barriered reverse `N2 + O -> NO + N` is measured in a frozen reservoir
    whose translational temperature is fixed while the internal
    (rotational/vibrational/electronic) modes are held at a different
    temperature (`in.reverse_noneq`, relaxation zeroed).  The ratio of the
    two-temperature reverse rate to the equilibrium rate at the same
    translational temperature is compared against an INDEPENDENT
    microcanonical integral of the reconstructed detailed-balance table over
    the two collision-energy distributions.  The rate swings ~14x with the
    internal temperature and matches the prediction to within statistics,
    directly verifying microscopic reversibility out of equilibrium.
15. **Nonlinear rotor**: an exchange pair `TRIA + ATB <-> DIA + ATO` (data in
    `nl.*`) whose reactant `TRIA` is a nonlinear triatomic (rotdof 3)
    exercises the nonlinear rotational partition function (`qrot ~ T^1.5`)
    and the `zcont = 3/2` continuum in the detailed-balance table.  Detailed
    balance holds to ~1% up to ~12 kK; a small (~4% at 20 kK) drift appears
    at very high temperature, beyond the table's built-in calibration
    self-check range - a documented limitation of the nonlinear-rotor path.
16. **Molecular third body under `vibrate smooth`**: check 13 with classical
    (continuum) vibration, so the third body's vibration is a flat measure
    variable folded into the 3-body density of states exactly as its rotation
    is; the recombination reproduces the classical-vibration volumetric
    `K_eq`.
17. **Sharply temperature-varying external `K_eq`**: a `keq_file` whose target
    differs from the statistical-mechanics `K_eq` by a factor that swings ~5x
    across the fit window stresses the Park fit of the residual
    `R(T)=K_eq_statmech/K_eq_ext`; the fit self-check must not warn and the
    measured reverse rate must reproduce the sharply varying target.
18. **Restart then continue**: a restart is written mid-run and read back, the
    `collide`/`react`/`fix` commands are re-issued (`in.reverse_restart1/2`),
    and the run continues.  The detailed-balance tables and `K_eq` fits are
    rebuilt deterministically at init and the per-particle electronic state is
    restored by `fix elecmode`, so the continued run holds detailed balance.
19. **Charge-exchange reverse** (`in.reverse_ce`): `MAp + MB <-> MA + MBp`
    moves charge between two atoms with no free electron, so it is an ordinary
    `EXCHANGE` whose reverse is derived by detailed balance (unlike ionization,
    whose reverse depends on the electron temperature and is rejected at init);
    the electronic ground-state degeneracies drive the `K_eq` prefactor.
20. **Multi-channel detailed balance** (`in.reverse_multi`): two coupled
    Zeldovich exchanges `N2+O<->NO+N` and `NO+O<->O2+N` share the reactants
    NO/O/N; both reverses are generated by `reverse auto` and each pair's tally
    ratio must independently equal its own analytic `K_eq` while competing for
    the shared species -- a test of reaction pairing/indexing under
    simultaneous active reverses.
21. **Multi-channel equilibrium benchmark** (`in.reverse_multieq`): the same two
    coupled pairs run in a closed reflective box seeded far from equilibrium;
    the five-species composition must relax to the joint chemical equilibrium,
    compared against an independent two-reaction-extent equilibrium solve (N and
    O nuclei conserved) rather than only internal self-consistency.
22. **Energy conservation across the reverse-disposal path**
    (`in.reverse_econsv`, `econsv.species`): a closed reacting box with the
    mass-consistent exchange `N2+O<->NO+N` (electronic off) dumps per-particle
    `ke+erot+evib` at the first and last step; the thermal-energy change must
    equal the net reaction extent times the reaction energy `coeff[4]` to
    machine precision (measured residual ~1e-8 of the thermal energy).  A
    disposal that leaked energy would fail by `O(reaction energy)` per event.
23. **Mode-resolved non-equilibrium** (`in.reverse_modeneq`): like check 14 but
    each internal mode (rotation, vibration, electronic) is driven hot
    *individually* while the other two stay at the translational temperature.
    The measured reverse-rate ratio to the fully-cold rate must match an
    independent microcanonical integral that puts only that mode hot, confirming
    each mode's energy contributes to the reverse rate per microscopic
    reversibility (measured dev <=1.5% for all three modes at 8->25 kK).
24. **Unbounded-probability diagnostic**: a barrierless reverse whose derived
    rate exceeds the collision rate saturates the per-collision probability
    (`react_prob>1`); the one-time diagnostic must fire to warn the user while
    the run still completes and tallies reactions (it is a warning, not a stop).

A probability capped at 1 limits the rate of that reaction.  Where a
*forward* reaction is capped, its reverse table is limited to match, so
the pair still relaxes to the right equilibrium (see `doc/src/react.rst`;
with `data/air.tce` this concerns NO + O and NO + N dissociation).  A
capped *reverse*, as in check 24, is not carried back into its forward
reaction.

Accurate reverse rates need accurate partition functions: supply a
`rotfile` with symmetry numbers and an `elecfile` with the low-lying
electronic levels for every species that participates in a reverse
reaction.  Omitting them silently biases `K_eq` (here by the factors
noted in check 5).
