# OpenMMForceFields PR #430 validation notebooks

Notebooks for validating [openmm/openmmforcefields#430](https://github.com/openmm/openmmforcefields/pull/430),
which adds multi-residue molecule support (proteins/peptides) to the SMIRNOFF
template generator, replaces NetworkX with OpenMM's template matcher, and
handles virtual sites spanning residue boundaries.

The goal is to verify that the OMMFFS template-generator pathway
(`SMIRNOFFTemplateGenerator` + `openmm.app.ForceField.createSystem()`) produces
the same energies and forces as a direct Interchange parameterization
(`OFFForceField.create_openmm_system()`).

## Notebooks

| Notebook | Purpose |
|----------|---------|
| `01_energy_comparison_protein.ipynb` | Per-component energy and force comparison for `test-ala-3.pdb` and `test-aa.pdb` at PDB and dynamics-perturbed geometries. Tolerances: 0.01 kcal/mol energy, 1e-5 relative force. |
| `02_virtual_sites_cross_residue.ipynb` | Adds custom virtual site parameters (TrivalentLonePair, MonovalentLonePair, BondCharge) and verifies vsite counts and energies match between pathways for protein-only and protein+ligand systems. |
| `03_performance_linear.ipynb` | Wall-time benchmark of Interchange vs OMMFFS for linear alkanes of increasing size (5-200 carbons). |
| `04_performance_branched_ring.ipynb` | Wall-time benchmark for MiniDrugBank molecules and synthetic fused-ring systems, comparing performance across molecule size and ring count. |

## OpenMM constraint bug

During this work we found a bug in OpenMM's `ForceField.createSystem()` where a
leaked loop variable caused constraints from one residue template to be
incorrectly applied to bonds in a different residue. Details and a minimal
reproducer are in `openmm_constraint_bug_issue.md`. The fix is
[openmm/openmm#5236](https://github.com/openmm/openmm/pull/5236).

## Environment

- Python 3.13, OpenMM 8.5.0-dev, openmmforcefields `whole-molecule` branch
- `openff_no_water-3.0.0-alpha0.offxml` force field
- `micromamba activate ommffs-dev`

## Authorship

These notebooks are mostly AI-generated (Claude Code) with human review and
modifications.
