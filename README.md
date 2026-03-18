# anything-goes

A repo for anything users/developers have done with OpenFF and want to share. 

We will approve basically any relevant PR. Your putting stuff here comes with no obligation to implement testing or provide support. Material here may be modified/overwritten by other people. Try to keep things under a few megabytes so we don't end up with an unmanageable git tree.

If you are in the org/have permissions and can merge without approval, please do! If not a maintainer will come around ASAP to approve your PR. You're free to merge any time GH isn't totally preventing you from doing so.

Provided entirely without support, endorsement, or warranty, either from the contributors or the Open Force Field Consortium/Initiative itself. This lack of warranty extends to the integrity of the summary table below.

## Directory Summary

| Directory / File | Description |
|---|---|
| [besmarts-fit/](besmarts-fit/) | Force field fitting with `besmarts` and `smee`/`descent`; SMARTS hierarchy labeling; parameter optimization against QCFractal datasets; hessian support. |
| [bespokefit_webserver/](bespokefit_webserver/) | Docker/Traefik deployment of `openff-bespokefit` as an HTTPS webservice; proof-of-concept for ASAP Discovery FEC pipelines on AWS. |
| [build_molecule_many_ways/](build_molecule_many_ways/) | Export single molecules, packed boxes, solvated systems, and mixtures to GROMACS/AMBER/LAMMPS using OpenFF Interchange and Packmol. |
| [combine-amber-components/](combine-amber-components/) | Load separate AMBER prmtop/inpcrd files and merge them into one `Interchange` object; OpenFF Interchange + OpenMM. |
| [demos_workshops/](demos_workshops/) | Index of links to external OpenFF workshop and demo materials (2022–2024 RDKit UGM, CCPBioSim, OMSF). No local code. |
| [diff-assigned-parameters/](diff-assigned-parameters/) | Compare SMIRKS parameter assignments between two OpenFF force fields using `ForceField.label_molecules`; highlights changed parameter IDs. |
| [evaluator-on-nrp-helm/](evaluator-on-nrp-helm/) | Run `openff-evaluator` on NRP Kubernetes via Dask and Helm; includes pod configs, persistent storage, and a client script. |
| [interchange-drivers/](interchange-drivers/) | Compare GROMACS and OpenMM potential energies for an octanol-hexane mixture via `Interchange` energy drivers; Sage force field. |
| [ionic-liquid-energies/](ionic-liquid-energies/) | Prepare a BMIM-Tf2N ionic liquid pair, evaluate energies across multiple MD engines, minimize with OpenMM; uses Interchange drivers. |
| [ionic_liquid_scripts/](ionic_liquid_scripts/) | Generate GROMACS inputs for ionic liquid (imidazolium/TFSI) and deep eutectic solvent systems at varied compositions; OpenFF Toolkit + Packmol. |
| [mdanalysis_pdb_loader/](mdanalysis_pdb_loader/) | Load PDB files with guessed bonds via MDAnalysis → RDKit → OpenFF Topology; workaround for PDB files lacking explicit connectivity. |
| [metal-param-example/](metal-param-example/) | Manually add library charges, bonds, angles, torsions, and vdW for a zinc organometallic using OpenFF Toolkit and Interchange. |
| [minimization-with-nagl.ipynb](minimization-with-nagl.ipynb) | Minimize multiple ligand conformers in a protein receptor; NAGL charges + Sage + ff14SB + OpenMM via `Interchange.combine()`. |
| [mix_openff_and_ommffs_components/](mix_openff_and_ommffs_components/) | Parameterize Mg²⁺ complexes by combining OpenMM amber/tip3p_HFE_multivalent for metal ions with Sage + ff14SB via Interchange. |
| [mnsol_assessment/](mnsol_assessment/) | Audit ~2400 MNSol solvation database entries for name/structure mismatches; XYZ→RDKit perception, graph isomorphism, OpenFF SMILES validation, TMOS. |
| [nagl-biopolymer/](nagl-biopolymer/) | Apply NAGL GNN charges to a solvated biopolymer/dye system with charge smearing for neutralization; `NAGLToolkitWrapper` + Interchange. |
| [nagl-ligand/](nagl-ligand/) | Minimal example: assign NAGL GNN AM1-BCC-style charges to a small molecule, create Interchange with Sage, export to MD engines. |
| [openff-3.0.0-alpha-notebook/](openff-3.0.0-alpha-notebook/) | Demo of OpenFF 3.0.0-alpha protein force field + OPC3 water on a solvated protein-ligand complex; short OpenMM simulation; MDTraj/NGLView. |
| [optimize-dimer/](optimize-dimer/) | Load a bromobenzene–water dimer from two SDF files, parameterize with Sage/Interchange, energy-minimize with OpenMM, compare coordinates. |
| [pdb_sdf_combined.ipynb](pdb_sdf_combined.ipynb) | Combine a protein PDB with an SDF ligand in OpenFF Toolkit; parameterize using AMBER ff14SB + GAFF via `GAFFTemplateGenerator`. |
| [pdbccdutils/](pdbccdutils/) | Load a complex PDB (nucleic acid + modified carbohydrate) via `pdbeccdutils` CLC reader, then simulate with OpenFF Toolkit. |
| [popc_two_charges/](popc_two_charges/) | Compare AM1-BCC (~12 min) vs. NAGL (~1 sec) charge assignment for POPC lipid parameterization; export GROMACS topology via Interchange. |
| [ptm_microviridin/](ptm_microviridin/) | Simulate post-translationally modified proteins (microviridin) using `openff-pablo` PDB loader, Sage + ff14SB + NAGL charges, and OpenMM. |
| [sage-230-conf-min/](sage-230-conf-min/) | Generate intentionally bad RDKit conformers for 50 drug fragments, then batch-minimize with Sage 2.3.0-rc2 via OpenFF Toolkit. |
| [simple-protein-ligand-nagl/](simple-protein-ligand-nagl/) | Minimal demo of NAGL GNN charges on a BRD4 protein-ligand complex; Sage + ff14SB; `NAGLToolkitWrapper` + `charge_from_molecules` via Interchange. |
| [system-generator-replacement/](system-generator-replacement/) | Side-by-side migration guide from `openmmforcefields.SystemGenerator` to native OpenFF Toolkit/Interchange for 4 common system types. |
