"""Shared utilities for the virtual-site export validation suite.

These helpers are used across notebooks 01-05.  Several are adapted from the
prior ``anything-goes/ommffs-pr430-tests/`` notebooks (energy decomposition,
parent-atom vsite matching) and extended with nonbonded/vsite introspection
needed to *diagnose* discrepancies rather than merely assert on totals.

Design rule (per Jeff): when a comparison fails, surface the discrepancy in
detail; never silently loosen a tolerance.
"""

from __future__ import annotations

import copy

import numpy as np
import openmm
import openmm.unit as unit

# ---------------------------------------------------------------------------
# Units and default tolerances
# ---------------------------------------------------------------------------
ENERGY_UNIT = unit.kilocalories_per_mole
FORCE_UNIT = unit.kilocalories_per_mole / unit.angstroms
LENGTH_UNIT = unit.nanometer

# "same force field, two export pathways" should agree to ~machine precision.
TIGHT_ENERGY_TOL = 1.0e-6 * ENERGY_UNIT      # kcal/mol, absolute-ish for small systems
TIGHT_FORCE_TOL = 1.0e-6                      # relative
TIGHT_POSITION_TOL = 1.0e-9                   # nm, vsite placement, same encoding
# cross-representation water (LocalCoordinatesSite vs ThreeParticleAverage/OutOfPlane)
CROSS_POSITION_TOL = 1.0e-6                   # nm
CROSS_ENERGY_TOL = 1.0e-4 * ENERGY_UNIT      # kcal/mol


# ---------------------------------------------------------------------------
# Energy / force decomposition
# ---------------------------------------------------------------------------
def compute_energy(system, positions, *, recompute_vsites=False, platform="Reference"):
    """Return per-force-group energies and forces for ``system`` at ``positions``.

    Each force is assigned a unique force group so components can be read back
    by class name.  ``recompute_vsites=True`` recomputes virtual-site positions
    from their parent atoms after setting positions -- use this whenever the
    supplied ``positions`` may not already contain consistent vsite coordinates
    (e.g. hand-built inputs, or positions borrowed from a differently-ordered
    system).  When positions already carry correct vsite coordinates (e.g. from
    an OpenMM integrator or ``Modeller.addExtraParticles``) it is a no-op.
    """
    system = copy.deepcopy(system)
    for index, force in enumerate(system.getForces()):
        force.setForceGroup(index)
    plat = openmm.Platform.getPlatformByName(platform)
    integrator = openmm.VerletIntegrator(0.001)
    context = openmm.Context(system, integrator, plat)
    if recompute_vsites:
        positions = _pad_positions(positions, system.getNumParticles())
    context.setPositions(positions)
    if recompute_vsites:
        context.computeVirtualSites()
    context.applyConstraints(integrator.getConstraintTolerance())

    def _group(i, **kw):
        return context.getState(groups=(1 << i), **kw)

    energy = {
        "total": context.getState(getEnergy=True).getPotentialEnergy(),
        "components": {
            system.getForce(i).__class__.__name__: _group(i, getEnergy=True).getPotentialEnergy()
            for i in range(system.getNumForces())
        },
    }
    forces = {
        "total": context.getState(getForces=True).getForces(asNumpy=True),
        "components": {
            system.getForce(i).__class__.__name__: _group(i, getForces=True).getForces(asNumpy=True)
            for i in range(system.getNumForces())
        },
    }
    del context, integrator
    return energy, forces


def _pad_positions(positions, n_particles):
    """Pad real-atom positions with zero placeholders for trailing vsites."""
    try:
        n_given = len(positions)
    except TypeError:
        n_given = positions.shape[0]
    if n_given >= n_particles:
        return positions
    arr = np.asarray(positions.value_in_unit(LENGTH_UNIT)) if hasattr(positions, "value_in_unit") \
        else np.asarray([[v.x, v.y, v.z] for v in positions])
    pad = np.zeros((n_particles - n_given, 3))
    return np.vstack([arr, pad]) * LENGTH_UNIT


def resolved_vsite_positions(system, positions, *, apply_constraints=False,
                             platform="Reference"):
    """Return the Nx3 numpy array (nm) of positions after OpenMM places vsites.

    Feeds ``positions`` (real-atom coordinates matter; vsite entries may be
    placeholders) into a context, optionally applies distance constraints (e.g.
    rigid-water geometry), calls ``computeVirtualSites``, and reads back all
    particle positions.  Used to check vsite geometry to machine precision.
    """
    system = copy.deepcopy(system)
    plat = openmm.Platform.getPlatformByName(platform)
    integrator = openmm.VerletIntegrator(0.001)
    context = openmm.Context(system, integrator, plat)
    positions = _pad_positions(positions, system.getNumParticles())
    context.setPositions(positions)
    if apply_constraints:
        context.applyConstraints(1e-10)
    context.computeVirtualSites()
    state = context.getState(getPositions=True)
    arr = state.getPositions(asNumpy=True).value_in_unit(LENGTH_UNIT)
    del context, integrator
    return np.asarray(arr)


# ---------------------------------------------------------------------------
# Particle bookkeeping: atoms vs virtual sites, parent-atom vsite matching
# ---------------------------------------------------------------------------
def count_vsites(system):
    return sum(system.isVirtualSite(i) for i in range(system.getNumParticles()))


def _vsite_parent_particles(system, vsite_idx):
    vs = system.getVirtualSite(vsite_idx)
    return tuple(vs.getParticle(j) for j in range(vs.getNumParticles()))


def real_atom_positions(system, positions):
    """Positions (Quantity, nm) of the non-vsite particles, in particle order."""
    positions = _pad_positions(positions, system.getNumParticles())
    if hasattr(positions, "value_in_unit"):
        arr = np.asarray(positions.value_in_unit(LENGTH_UNIT))
    else:
        arr = np.asarray([[v.x, v.y, v.z] for v in positions])
    idx = [i for i in range(system.getNumParticles()) if not system.isVirtualSite(i)]
    return arr[idx] * LENGTH_UNIT


def build_atom_maps(system):
    """Map particle indices to atom ordinals (vsites -> None)."""
    particle_to_atom, atom_to_particle = {}, {}
    atom_ord = 0
    for i in range(system.getNumParticles()):
        if system.isVirtualSite(i):
            particle_to_atom[i] = None
        else:
            particle_to_atom[i] = atom_ord
            atom_to_particle[atom_ord] = i
            atom_ord += 1
    return particle_to_atom, atom_to_particle


def get_permutation_indices(system_1, system_2):
    """Permutation mapping particles system_1 <-> system_2.

    Real atoms are matched by sequential order (both stacks keep molecules
    contiguous).  Virtual sites are matched by the *ordinals of their parent
    atoms*, not by sequential vsite order -- Interchange puts all vsites at the
    end of the system while OpenMM/openmmforcefields interleave them per
    residue.  Returns ``(perm_12, perm_21)`` where ``perm_12[i]`` is the
    system_2 index corresponding to system_1 index ``i``.
    """
    atoms_1 = [i for i in range(system_1.getNumParticles()) if not system_1.isVirtualSite(i)]
    sites_1 = [i for i in range(system_1.getNumParticles()) if system_1.isVirtualSite(i)]
    atoms_2 = [i for i in range(system_2.getNumParticles()) if not system_2.isVirtualSite(i)]
    sites_2 = [i for i in range(system_2.getNumParticles()) if system_2.isVirtualSite(i)]

    assert len(atoms_1) == len(atoms_2), f"Atom count mismatch: {len(atoms_1)} vs {len(atoms_2)}"
    assert len(sites_1) == len(sites_2), f"Vsite count mismatch: {len(sites_1)} vs {len(sites_2)}"

    atom_to_ord_1 = {p: o for o, p in enumerate(atoms_1)}
    atom_to_ord_2 = {p: o for o, p in enumerate(atoms_2)}

    def parent_ords(system, vi, a2o):
        return tuple(a2o[p] for p in _vsite_parent_particles(system, vi))

    sites2_by_parents = {}
    for si in sites_2:
        sites2_by_parents.setdefault(parent_ords(system_2, si, atom_to_ord_2), []).append(si)

    site_pairs = []
    for si in sites_1:
        key = parent_ords(system_1, si, atom_to_ord_1)
        matched = sites2_by_parents.get(key, [])
        assert matched, f"No matching vsite in system_2 for vsite {si} (parents {key})"
        site_pairs.append((si, matched.pop(0)))

    n = system_1.getNumParticles()
    perm_12, perm_21 = [-1] * n, [-1] * n
    for a1, a2 in zip(atoms_1, atoms_2):
        perm_12[a1], perm_21[a2] = a2, a1
    for s1, s2 in site_pairs:
        perm_12[s1], perm_21[s2] = s2, s1
    return perm_12, perm_21


# ---------------------------------------------------------------------------
# Virtual-site geometry introspection
# ---------------------------------------------------------------------------
def vsite_details(system, vsite_idx):
    """Type + geometry parameters for a vsite particle."""
    vs = system.getVirtualSite(vsite_idx)
    info = {"type": vs.__class__.__name__,
            "parents": list(_vsite_parent_particles(system, vsite_idx))}
    if isinstance(vs, openmm.LocalCoordinatesSite):
        lp = vs.getLocalPosition()
        info.update(
            originWeights=list(vs.getOriginWeights()),
            xWeights=list(vs.getXWeights()),
            yWeights=list(vs.getYWeights()),
            localPosition=[lp.x, lp.y, lp.z],
        )
    elif isinstance(vs, openmm.TwoParticleAverageSite):
        info["weights"] = [vs.getWeight(0), vs.getWeight(1)]
    elif isinstance(vs, openmm.ThreeParticleAverageSite):
        info["weights"] = [vs.getWeight(0), vs.getWeight(1), vs.getWeight(2)]
    elif isinstance(vs, openmm.OutOfPlaneSite):
        info["weights"] = [vs.getWeight12(), vs.getWeight13(), vs.getWeightCross()]
    return info


# ---------------------------------------------------------------------------
# Nonbonded introspection (per-particle params + exceptions)
# ---------------------------------------------------------------------------
def _nonbonded_force(system):
    forces = [f for f in system.getForces() if isinstance(f, openmm.NonbondedForce)]
    if len(forces) != 1:
        raise ValueError(
            f"Expected exactly one NonbondedForce, found {len(forces)}. "
            "This system likely used combine_nonbonded_forces=False; use "
            "nonbonded_particle_params / split-force helpers instead."
        )
    return forces[0]


def particle_nb_params(system, idx):
    """(charge_e, sigma_A, epsilon_kcal) for a particle from the NonbondedForce."""
    q, sig, eps = _nonbonded_force(system).getParticleParameters(idx)
    return {
        "charge": q.value_in_unit(unit.elementary_charge),
        "sigma": sig.value_in_unit(unit.angstrom),
        "epsilon": eps.value_in_unit(unit.kilocalories_per_mole),
    }


def exceptions_table(system):
    """All NonbondedForce exceptions keyed by sorted particle-index pair.

    Value: dict with chargeProd (e^2), sigma (A), epsilon (kcal/mol) and a
    boolean ``is_exclusion`` (all three ~0).
    """
    nb = _nonbonded_force(system)
    table = {}
    for ei in range(nb.getNumExceptions()):
        p1, p2, q, sig, eps = nb.getExceptionParameters(ei)
        key = (min(p1, p2), max(p1, p2))
        qv = q.value_in_unit(unit.elementary_charge ** 2)
        sv = sig.value_in_unit(unit.angstrom)
        ev = eps.value_in_unit(unit.kilocalories_per_mole)
        table[key] = {
            "chargeProd": qv, "sigma": sv, "epsilon": ev,
            "is_exclusion": abs(qv) < 1e-12 and abs(ev) < 1e-12,
        }
    return table


def total_charge(system):
    nb = _nonbonded_force(system)
    return sum(
        nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)
        for i in range(system.getNumParticles())
    )


# ---------------------------------------------------------------------------
# First-principles (graph-derived) exception checks
# ---------------------------------------------------------------------------
def graph_distances(n_atoms, bonds):
    """All-pairs shortest-path distance (in bonds) via BFS. ``bonds`` is a list
    of ``(i, j)`` real-atom index pairs. Returns an ``n x n`` list of ints
    (unreachable pairs are a large sentinel)."""
    adj = {i: set() for i in range(n_atoms)}
    for a, b in bonds:
        adj[a].add(b)
        adj[b].add(a)
    INF = 10 ** 9
    dist = [[INF] * n_atoms for _ in range(n_atoms)]
    for s in range(n_atoms):
        dist[s][s] = 0
        frontier, seen, d = [s], {s}, 0
        while frontier:
            nxt = []
            for u in frontier:
                for w in adj[u]:
                    if w not in seen:
                        seen.add(w)
                        dist[s][w] = d + 1
                        nxt.append(w)
            frontier, d = nxt, d + 1
    return dist


def exception_kinds(system):
    """``{(i,j): "excl"|"scaled"}`` for every NonbondedForce exception (sorted
    key). Pairs absent from the table are implicitly full-strength."""
    out = {}
    for (a, b), e in exceptions_table(system).items():
        out[(min(a, b), max(a, b))] = "excl" if e["is_exclusion"] else "scaled"
    return out


def vsite_parent0_and_frame(system, vsite_idx):
    """(orientation-atom-0 particle index, set of all frame particle indices).

    Interchange registers each vsite under ``orientations[0]`` and OpenMM's
    app layer shares its exclusions with the first frame particle, so a vsite
    inherits atom-0's exception relationships and is additionally excluded from
    every atom in its own frame."""
    vs = system.getVirtualSite(vsite_idx)
    parents = [vs.getParticle(j) for j in range(vs.getNumParticles())]
    return parents[0], set(parents)


def full_exception_values(system):
    """``{parent-ordinal key: (chargeProd, sigma, epsilon)}`` for every
    exception, keyed so tables from two differently-ordered systems (Interchange
    puts vsites last; openmmforcefields interleaves) can be compared directly.
    Atom keys are ``("A", ordinal)``; vsite keys are ``("V", sorted parent
    ordinals)``."""
    p2a, _ = build_atom_maps(system)

    def label(i):
        if system.isVirtualSite(i):
            vs = system.getVirtualSite(i)
            pars = tuple(sorted(p2a[vs.getParticle(j)] for j in range(vs.getNumParticles())))
            return ("V", pars)
        return ("A", p2a[i])

    table = {}
    for (a, b), e in exceptions_table(system).items():
        key = tuple(sorted([label(a), label(b)], key=str))
        table[key] = (round(e["chargeProd"], 8), round(e["sigma"], 6), round(e["epsilon"], 8))
    return table


def compare_exception_tables(system_a, system_b, *, tol=1e-6):
    """Diff two systems' full exception tables (values, parent-matched). Returns
    a list of ``(key, values_a, values_b)`` for pairs that are missing on one
    side or whose (chargeProd, sigma, epsilon) differ by more than ``tol``."""
    ta, tb = full_exception_values(system_a), full_exception_values(system_b)
    diffs = []
    for key in sorted(set(ta) | set(tb), key=str):
        va, vb = ta.get(key), tb.get(key)
        if va is None or vb is None or any(abs(x - y) > tol for x, y in zip(va, vb)):
            diffs.append((key, va, vb))
    return diffs


# ---------------------------------------------------------------------------
# High-level comparison
# ---------------------------------------------------------------------------
def _norm_sq(x):
    return (x * x).sum() / x.shape[0]


def relative_force_deviation(f1, f2):
    a = f1.value_in_unit(FORCE_UNIT)
    b = f2.value_in_unit(FORCE_UNIT)
    scale = _norm_sq(a) + _norm_sq(b)
    return float(np.sqrt(_norm_sq(a - b) / scale)) if scale > 0 else 0.0


def compare_energies(name, positions, template_system, reference_system,
                     *, energy_tol=TIGHT_ENERGY_TOL, force_tol=TIGHT_FORCE_TOL,
                     recompute_vsites=False, verbose=True):
    """Compare per-component energies + forces between two systems.

    ``positions`` are in ``template_system`` ordering; they are permuted into
    ``reference_system`` ordering via parent-atom vsite matching.  Raises
    AssertionError (with a detailed printout) if totals disagree beyond tol.
    """
    forces_perm, pos_perm = get_permutation_indices(template_system, reference_system)

    t_e, t_f = compute_energy(template_system, positions, recompute_vsites=recompute_vsites)
    ref_positions = [positions[i] for i in pos_perm]
    r_e, r_f = compute_energy(reference_system, ref_positions, recompute_vsites=recompute_vsites)

    r_f["total"] = np.array(
        [r_f["total"][i].value_in_unit(FORCE_UNIT) for i in forces_perm]
    ) * FORCE_UNIT

    delta = t_e["total"] - r_e["total"]
    rel_dev = relative_force_deviation(t_f["total"], r_f["total"])

    if verbose:
        comps = set(r_e["components"]) & set(t_e["components"])
        print(f"\n{'Component':<28} {'Template':>18} {'Reference':>18} {'Delta':>14}")
        print("-" * 80)
        for key in sorted(comps):
            te = t_e["components"][key].value_in_unit(ENERGY_UNIT)
            re = r_e["components"][key].value_in_unit(ENERGY_UNIT)
            print(f"{key:<28} {te:18.8f} {re:18.8f} {te - re:14.6e}")
        tt = t_e["total"].value_in_unit(ENERGY_UNIT)
        rt = r_e["total"].value_in_unit(ENERGY_UNIT)
        print("-" * 80)
        print(f"{'TOTAL':<28} {tt:18.8f} {rt:18.8f} {tt - rt:14.6e}")
        print(f"Relative force deviation: {rel_dev:.3e}")

    assert abs(delta) <= energy_tol, f"[{name}] energy deviation {delta} > {energy_tol}"
    assert rel_dev <= force_tol, f"[{name}] force deviation {rel_dev:.3e} > {force_tol}"
    if verbose:
        print(f"PASS: {name}")
    return delta, rel_dev


def diff_systems(sys_a, sys_b, name_a="A", name_b="B", *, tol=1e-6):
    """Diagnostic dump of per-atom/vsite nonbonded + vsite-geometry differences.

    Matches atoms by ordinal and vsites by parent-atom ordinals.  Prints only
    differences.  Returns the number of differences found (0 == identical).
    """
    p2a_a, a2p_a = build_atom_maps(sys_a)
    p2a_b, a2p_b = build_atom_maps(sys_b)
    ndiff = 0

    # real atoms
    for o in range(len(a2p_a)):
        na = particle_nb_params(sys_a, a2p_a[o])
        nb = particle_nb_params(sys_b, a2p_b[o])
        for k in ("charge", "sigma", "epsilon"):
            if abs(na[k] - nb[k]) > tol:
                ndiff += 1
                print(f"  atom ord {o}: {k} {name_a}={na[k]:.6f} {name_b}={nb[k]:.6f}")

    # vsites matched by parent ordinals
    def by_parents(system, p2a):
        d = {}
        for i in range(system.getNumParticles()):
            if system.isVirtualSite(i):
                key = tuple(p2a[p] for p in _vsite_parent_particles(system, i))
                d.setdefault(key, []).append(i)
        return d

    va, vb = by_parents(sys_a, p2a_a), by_parents(sys_b, p2a_b)
    for key in sorted(set(va) | set(vb), key=str):
        la, lb = va.get(key, []), vb.get(key, [])
        if len(la) != len(lb):
            ndiff += 1
            print(f"  vsite(parents {key}): count {name_a}={len(la)} {name_b}={len(lb)}")
            continue
        for ia, ib in zip(la, lb):
            na, nb = particle_nb_params(sys_a, ia), particle_nb_params(sys_b, ib)
            for k in ("charge", "sigma", "epsilon"):
                if abs(na[k] - nb[k]) > tol:
                    ndiff += 1
                    print(f"  vsite(parents {key}): {k} {name_a}={na[k]:.6f} {name_b}={nb[k]:.6f}")
    if ndiff == 0:
        print(f"  no nonbonded differences between {name_a} and {name_b}")
    return ndiff


# ---------------------------------------------------------------------------
# Dynamics
# ---------------------------------------------------------------------------
def propagate_dynamics(positions, system, nsteps=100, *, minimize=True,
                       temperature=300 * unit.kelvin, timestep=1.0 * unit.femtoseconds):
    integrator = openmm.LangevinMiddleIntegrator(
        temperature, 1.0 / unit.picoseconds, timestep
    )
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName("Reference"))
    context.setPositions(positions)
    if minimize:
        openmm.LocalEnergyMinimizer.minimize(context)
    integrator.step(nsteps)
    pos = context.getState(getPositions=True).getPositions()
    del context, integrator
    return pos
