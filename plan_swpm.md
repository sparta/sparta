# SWPM Future Work

## 1. Gradual / Continuous Reduction

### Motivation

The current scheme is all-or-nothing: a group must reach `Ncmax` before any
reduction fires.  In practice this means cells near the threshold see large
bursts of particle deletion followed by quiet periods, introducing temporal
noise in macroscopic quantities.  A gradual scheme applies light reduction
continuously so the particle count stays closer to a target without large
discrete jumps.

### Idea

Replace the hard `np > Ncmax` trigger with a soft probability that grows with
group size and shrinks as the count approaches `Ngmin`:

```
P_reduce(np) = clamp((np - N_target) / (Ncmax - N_target), 0, 1)^alpha
```

- `N_target` is the desired equilibrium count (between `Ngmin` and `Ncmax`).
- `alpha > 1` makes the ramp concave (slow to trigger near target, aggressive
  near `Ncmax`).  `alpha = 1` is linear.
- When `np <= N_target`, `P_reduce = 0`: no reduction fires.
- When `np = Ncmax`, `P_reduce = 1`: reduction always fires (recovers current
  behavior as a limit).

At each timestep, for each group that exceeds `N_target`, draw a uniform
random number and apply reduction only if it is less than `P_reduce(np)`.
This means on average only a fraction of over-threshold groups are reduced per
step, spreading the work in time.

### Strength scaling

Within a triggered reduction, the number of survivors can itself be
size-dependent rather than fixed at 2 (ENERGY/HEAT) or 2*nK (STRESS):

- Below some `N_soft`, emit the usual fixed number of survivors.
- Above `N_soft`, emit `floor(np * f_survive(np))` survivors, where
  `f_survive` decreases toward `2/np` as `np -> Ncmax`.

This lets heavily over-loaded cells shed particles aggressively while cells
just above target barely change shape.

### Implementation notes

- Needs two new `collide_modify` params: `N_target` and `alpha`
  (or `ramp_exponent`).
- `group_bt` already recurses; the probability check can sit in `group()`
  before dispatching to `reduce_*`.
- Conservation must still hold per triggered reduction — the gradual part is
  only in *whether* to reduce, not *how*.

---

## 2. Multi-species Support

### Motivation

All current SWPM logic (split, group, reduce) assumes a single species: the
stochastic weight is a scalar per particle and the grouping uses a single
velocity-space covariance.  Real flows (e.g. air = N2 + O2 + O + N) have
species-specific mass, and mixing species in one group biases the covariance
and misidentifies the thermal spread.

### Key issues

**Collision weight**  
`split()` sets `Gwtf = min(isw, jsw)`, which is species-agnostic.  For
cross-species collisions the weight should account for the species-weight ratio
(`specwt[isp] / specwt[jsp]`); otherwise rare-species particles are
preferentially split away.

**Grouping**  
`group_bt` builds a covariance matrix from all particles in the cell.  With
multiple species the covariance conflates thermal motion and species drift:
light atoms (O) have larger thermal velocities than heavy molecules (N2) at the
same temperature, so the eigenvectors of the mixed covariance do not align with
either species' true distribution.  The fix is to group per-species or to
mass-weight the covariance:

```
pij_ms[d1][d2] = sum_k  m_k * sw_k * (v_k[d1] - V[d1]) * (v_k[d2] - V[d2])
```

where `m_k` is the particle mass and `V` is the mass-weighted CoM velocity.

**Reduction conservation targets**  
`reduce_energy` / `reduce_heat` / `reduce_stress` conserve mass, momentum, and
energy computed from a uniform mass assumption.  With mixed species the
conserved quantities must be computed per species and the survivors must be
chosen from the correct species pool.

### Proposed approach

1. **Per-species grouping**: after `group_bt` partitions by velocity space,
   further split each leaf group by species before calling `reduce_*`.  This
   guarantees survivors are the same species as the particles they replace.

2. **Mass-weighted covariance in `group_bt`**: replace `pmsw` (momentum-weight)
   with `pmsw * mass[isp]` so the eigenvectors reflect actual thermal structure.

3. **Species-aware split**: in `split()`, check species equality; for
   cross-species pairs use the species weight ratio to set `Gwtf` correctly.

4. **Conserved-quantity accounting**: pass per-species mass, momentum, energy,
   Erot, Evib vectors to `reduce_*` and pick survivors from species-matched
   sub-lists.

### Open questions

- Should the grouping tree be built globally (then split by species at leaves)
  or independently per species?  The former reuses the tree structure; the
  latter is simpler but may produce very small per-species groups in trace
  species cells.
- How to handle a species with only 1 particle in a leaf group — skip reduction
  or merge into the nearest species by mass?
