# Virtual-site export validation suite

Validates that **OpenFF Interchange exports SMIRNOFF virtual sites to OpenMM
correctly**, both directly (`ForceField.create_openmm_system`) and through
`openmmforcefields`' `SMIRNOFFTemplateGenerator`.

The core representational fact under test: Interchange exports *every* SMIRNOFF
virtual site as an OpenMM `LocalCoordinatesSite`, whereas the canonical
OpenMM/AMBER water models use `OutOfPlaneSite` / `ThreeParticleAverageSite`. This
suite confirms the two agree on geometry, parameters, energies, and dynamics, and
regression-locks a cluster of virtual-site bugs found in early 2026.

## Notebooks

| Notebook | Theme | What it checks |
|----------|-------|----------------|
| `00_environment_and_provenance` | — | pins/records versions, gates minimums, writes `provenance.json` |
| `01_interchange_vs_openmmforcefields` | A | 4 vsite types × 4 charge models (library/user/AM1-BCC/NAGL); molecule reordering & intra-molecule remap; duplicate & mixed topologies; **full exception-table equality**; **non-default distinct 1-4 scaling**; split-force layout **(NoCutoff and PME)**; **robustness guards** (out-of-model settings must raise) |
| `02_openff_vs_external_water_models` | B | TIP5P / TIP4P-Ew / OPC: placement, charges, LJ, exclusions, PME box energy; random-orientation equivariance |
| `03_handcalc_exact_and_exclusions` | C | closed-form geometry for all 4 types incl. **general out-of-plane MonovalentLonePair**; charge-increment bookkeeping; 1-4 scaling incl. non-default; **exclusion set vs the bond graph** (single & multi-molecule); **end-to-end hand-summed Coulomb + LJ energy** with scaled vsite 1-4 pairs |
| `04_dynamics_stability` | D | OPC/TIP5P water boxes and virtual sites in *both* solute and solvent: stability, temperature, density; IC vs OMMFFS box energy |
| `05_findings_summary` | — | consolidated pass/fail table and findings |
| `07_diverse_box_both_pathways` | D | an **adversarial pathway-equivalence** box: all 12 diverse vsite solutes in **mostly OPC water**, **two copies each placed noncontiguously** with **one copy atom-reversed** (`Molecule.remap`), and vsite definitions written **inline with every charge-increment/LJ value unique & nonzero** (incl. C–F, and two solutes carrying two same-type `BondCharge` vsites with different params). Looped over **three charge methods** — `NAGLCharges`, `ToolkitAM1BCC`, and **user charges** (`charge_from_molecules`, which also stresses charge mapping onto the reversed copies). Asserts, per method: reversed copy ≡ in-order copy, box net-neutral, **Interchange vs `openmmforcefields` single-point energies agree to machine precision** (~1e-7), and stable NPT (HMR/4 fs) both ways; writes a PDB + DCD per (method × pathway) and reports T & density (± std dev). Few functions, no `_helpers` |
| `06_sfe_representative_vsite_dynamics` | D | **real freesolv/mnsol solutes** covering all 4 vsite types (carbonyl, aromatic-N, halogen σ-hole, amine-N) **plus two multi-vsite solutes** (each carrying two types at once): expected-parameter checks **through both export pathways** + **~1 ns HMR/4 fs NPT MD** (full 0.9 nm PME cutoff, HMR via each toolkit's native argument) in OPC water with vsites on solute **and** solvent (single-solute, multi-vsite solute, mixed multi-type box) across **three distinct charge mechanisms** — `NAGLCharges`/`ToolkitAM1BCC` in the force field and **user charges** (`charge_from_molecules`) — with the charges *verified different* (+ a scaled-charge positive control); parameters as expected, Interchange vs `openmmforcefields` energies agree to machine precision, systems don't explode; writes inspectable trajectories (`traj_06/`) with an `nglview` `view()` helper |

Helper modules: `_helpers.py` (energy decomposition, parent-atom vsite matching,
nonbonded/geometry introspection, **bond-graph exclusion derivation & exception-
table diffing**), `pathways.py` (the two export routes + a deterministic vsite
force-field builder).

Notebook 06 also reads `vsites-regression.offxml` — OpenFF's maintained all-types
test virtual-site block, vendored verbatim from
`interchange-regression-testing/openmm-system-parity/force-fields/virtual-sites.offxml`
(`BondCharge` on C–Cl/C–Br, `MonovalentLonePair` on carbonyls, `TrivalentLonePair`
on amine N, `DivalentLonePair` on aromatic-ring N). Notebook 06 grafts that block
onto Sage 2.3 + OPC and adds a C–F `BondCharge`. Its solute SMILES are drawn from
the freesolv/mnsol solvation-free-energy sets in
`github.com/openforcefield/ash-sage-rc2/tree/main/04_benchmark/sfes/input-data`,
plus two constructed **multi-vsite** solutes (4'-chloroacetophenone —
`MonovalentLonePair` + `BondCharge`; 3-bromopyridine — `DivalentLonePair` +
`BondCharge`) that exercise two virtual-site types on a single molecule.
The virtual-site numeric parameters are deliberately placeholders — this notebook
checks that the export/dynamics *machinery* is correct and stable on real
benchmark molecules, not that the parameters are physically optimal (those will
come from the fitted virtual-site force field).

Each box is packed just large enough to use the full 0.9 nm PME cutoff with
margin for the barostat, and propagated ~1 ns at a 4 fs timestep with hydrogen-mass
repartitioning requested through each toolkit's **native** argument (Interchange's
`to_openmm_system(hydrogen_mass=…)` and OpenMM's `createSystem(hydrogenMass=…)`).
Section A checks parameter assignment through
**both** export pathways; section B runs the `openmmforcefields` pathway (including
the two multi-vsite solutes); sections
C and D build **both** pathways from identical coordinates, assert their potential
energies agree to machine precision, and then run independent NPT dynamics through
each. Section D exercises three **distinct** charge mechanisms on one shared box —
`NAGLCharges` in the force field, `ToolkitAM1BCC` in the force field, and **user
charges** injected via `charge_from_molecules` (a separate code path) — and *proves*
they differ: it tabulates the per-atom solute charges, asserts `NAGLCharges` ≠
`ToolkitAM1BCC` (they differ by only ~0.006 e because NAGL is a GNN trained to
reproduce AM1-BCC, which is why their box total energies are nearly identical), and
runs a scaled-charge positive control showing the single-point energy moves ~36×
more when the charges genuinely change. Every simulation streams a topology + 100
frames (including virtual-site particles) to `traj_06/`; after running the
notebook, `view("<name>")` renders any of them with `nglview`. Running notebook 06 also installs `mdtraj` and `nglview`
(now pinned in `environment.yml`). A fast dry run is available by setting
`NB06_PROD_NS`, `NB06_EQUIL_NS`, `NB06_N_WATER`, and `NB06_N_FRAMES`.

### What the two pathways in notebook 01 do and don't prove

`openmmforcefields`' `SMIRNOFFTemplateGenerator` builds its `System` by calling
`ForceField.create_openmm_system`, which on this toolkit **is Interchange**; it
then serializes to an `ffxml` and lets `openmm.app.ForceField` rebuild it. The
two routes therefore **share Interchange's virtual-site math** — a bug inside it
would appear in both and cancel. Notebook 01 validates the `openmmforcefields`-
specific `ffxml` round-trip (parameters, charges, and the full exception table
regenerated from `coulomb14scale`/`lj14scale` and single-parent vsite sharing).
Interchange's virtual-site math is validated against **independent** references
in notebook 02 (native OpenMM water models) and notebook 03 (closed-form geometry,
bond-graph exclusion sets, and hand-summed nonbonded energies).

## Environment

Pinned, released versions (see `environment.yml`):

```
micromamba create -f environment.yml
micromamba activate vsite-validation
```

Key pins: `openmm=8.5.2`, `openmmforcefields=0.16.0` (first release with SMIRNOFF
virtual-site + multi-residue support), `openforcefields=2026.1.0`,
`openff-nagl-models=2025.9.0` (bundles the production `openff-gnn-am1bcc-1.0.0`
model used by Sage 2.3), `ambertools`, and `packmol` (solvation in notebook 04).

**Interchange is installed from `main`**, because the fix for #1498 was merged
after the 0.5.3 release:

```
pip install --no-deps 'git+https://github.com/openforcefield/openff-interchange.git@main'
```

Notebook 00 asserts `openff-interchange > 0.5.3` so the environment cannot
silently run without the #1498 fix.

## Running

```
micromamba run -n vsite-validation python build_notebooks.py --all   # .py -> .ipynb
for nb in 0*.ipynb; do
  micromamba run -n vsite-validation jupyter nbconvert --to notebook --execute \
    --inplace --ExecutePreprocessor.timeout=1800 "$nb"
done
```

Every numerical check is an `assert`, so a clean execution is a pass. The `.py`
files are the source of truth (py:percent format); `build_notebooks.py` converts
them to notebooks.

## Virtual-site bugs regression-locked

All surfaced in early 2026 while cross-checking Interchange against
`openmmforcefields`, fixed in Interchange 0.5.1 unless noted.

| Issue | Summary | Fixed in | Exercised by |
|-------|---------|----------|--------------|
| interchange #1436 | non-default 1-4 scaling not applied to vsite exceptions | 0.5.1 | `03`, `01` |
| interchange #1438 | `DivalentLonePair` geometry wrong | 0.5.1 | `03` |
| interchange #939  | `MonovalentLonePair` position off | (tracked) | `03` |
| interchange #1498 | vsite exceptions zeroed with `combine_nonbonded_forces=False` | after 0.5.3 (`main`) | `01` (fix verified) |

## Findings

**All 65 checks pass** on the certified environment. Interchange's SMIRNOFF
virtual-site export — geometry, charges, exclusion sets, 1-4 scaling, and total
nonbonded energies — matches independent closed-form and hand-summed references
to machine precision, both directly and through `openmmforcefields`, including
the split-force and PME layouts.

Notebook 01 verifies that the #1498 fix (merged after the 0.5.3 release) is
correct and complete: with `combine_nonbonded_forces=False` the split- and
combined-force total energies (and forces) now agree for every virtual-site type,
in a periodic PME box as well as NoCutoff, and a no-vsite control confirms the
split-force machinery is sound. On the unpatched 0.5.3 this same case was a
6.14 kcal/mol error for a `BondCharge` system with scaled 1-4 vsite exceptions;
the default combined-force path used by essentially all consumers was always
correct.

### No silent-danger input values

Force-field settings that would violate OpenMM's hardwired exception model do not
silently produce wrong energies — they raise:

| Setting | Rejected by | When |
|---------|-------------|------|
| `scale12` ≠ 0, `scale13` ≠ 0, `scale15` ≠ 1 (electrostatics or vdW) | OpenFF Toolkit | force-field parse |
| `exclusion_policy` ≠ `"parents"` | Interchange (and `openmmforcefields`) | export |
| `scale14` (any value) | — accepted, applied consistently by both routes | — |

Both Interchange and OpenMM's app layer inherit a virtual site's nonbonded
exclusions from a **single** parent (its first orientation atom / `excludeWith`),
cloning that atom's exception set; the two agree, and notebook 03 confirms the
resulting set matches a first-principles bond-graph derivation.

## Deferred follow-up

Copying `LocalCoordinatesSite` virtual sites into a **hybrid-topology** alchemical
`System` (for relative free-energy calculations) is a separate downstream concern
(historically openfe #394/#395, written before Interchange emitted
`LocalCoordinatesSite`). A follow-up notebook should build the hybrid `System` via
the alchemical factory and confirm virtual sites survive with correct parents,
exceptions, and end-state energies. The export foundations it relies on are
validated here.
