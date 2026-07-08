"""The two export pathways under test, plus deterministic vsite force-field builders.

Pathway 1 (reference): OpenFF Toolkit -> Interchange -> OpenMM System, via
``ForceField.create_openmm_system``.

Pathway 2: openmmforcefields ``SMIRNOFFTemplateGenerator`` registered on an
``openmm.app.ForceField``, then ``createSystem`` on an OpenMM topology whose
extra particles were added by ``Modeller.addExtraParticles``.

Both consume the *same* SMIRNOFF force field, so a correct implementation must
produce identical parameters/energies (modulo particle ordering, which
``_helpers.get_permutation_indices`` reconciles).
"""

from __future__ import annotations

import io

import openmm
from openff.toolkit import ForceField as OFFForceField
from openmm.app import ForceField as OMMForceField
from openmm.app import Modeller, NoCutoff

from openmmforcefields.generators import SMIRNOFFTemplateGenerator

# SMIRNOFF 0.4 tag headers copied verbatim from a shipped modern force field
# (openff-2.2.0.offxml) so custom test force fields parse under toolkit 0.18.
_BONDS = ('<Bonds version="0.4" potential="harmonic" '
          'fractional_bondorder_method="AM1-Wiberg" '
          'fractional_bondorder_interpolation="linear">'
          '<Bond smirks="[*:1]~[*:2]" length="1.5*angstrom" '
          'k="600.0*angstrom**-2*mole**-1*kilocalorie"/></Bonds>')
_ANGLES = ('<Angles version="0.3" potential="harmonic">'
           '<Angle smirks="[*:1]~[*:2]~[*:3]" angle="109.5*degree" '
           'k="120.0*mole**-1*radian**-2*kilocalorie"/></Angles>')
_TORS = ('<ProperTorsions version="0.4" '
         'potential="k*(1+cos(periodicity*theta-phase))" default_idivf="auto">'
         '<Proper smirks="[*:1]~[*:2]~[*:3]~[*:4]" periodicity1="3" '
         'phase1="0*degree" k1="0.0*mole**-1*kilocalorie" idivf1="1.0"/>'
         '</ProperTorsions>')


def library_charge_block(library_charges):
    """Render a ``<LibraryCharges>`` block from ``(smirks, charge)`` pairs."""
    lib = "".join(
        f'<LibraryCharge smirks="{s}" charge1="{q}*elementary_charge"/>'
        for s, q in library_charges
    )
    return f'<LibraryCharges version="0.3">{lib}</LibraryCharges>'


def build_offxml(vdw_atoms, library_charges=None, vsites=(), *, charge_block=None,
                 elec_scale14=0.833333, vdw_scale14=0.5,
                 elec_scale12=0.0, elec_scale13=0.0, elec_scale15=1.0,
                 vdw_scale12=0.0, vdw_scale13=0.0, vdw_scale15=1.0,
                 exclusion_policy="parents"):
    """Assemble a deterministic vsite force field as a SMIRNOFF string.

    ``vdw_atoms``: list of ``(smirks, sigma_A, epsilon_kcal)``.  The base charge
    handler is either ``charge_block`` (a raw ``<LibraryCharges>`` /
    ``<ToolkitAM1BCC>`` / ``<ChargeIncrementModel>`` XML string) or, if that is
    ``None``, a ``<LibraryCharges>`` block generated from ``library_charges``.
    ``vsites``: list of raw ``<VirtualSite .../>`` XML strings.
    ``elec_scale14``/``vdw_scale14`` drive non-default 1-4 scaling (#1436).
    The ``*_scale12``/``*_scale13``/``*_scale15`` and ``exclusion_policy``
    knobs exist only so the robustness section (notebook 01) can show that
    non-default adjacent scalings and unsupported vsite policies are *rejected*
    (never silently mis-handled); every real force field leaves them at the
    defaults printed here.
    """
    vdw = "".join(
        f'<Atom smirks="{s}" epsilon="{e}*kilocalorie/mole" sigma="{sig}*angstrom"/>'
        for s, sig, e in vdw_atoms
    )
    if charge_block is None:
        charge_block = library_charge_block(library_charges or [])
    vs = "".join(vsites)
    return f"""<?xml version="1.0" encoding="ASCII"?>
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
{_BONDS}
{_ANGLES}
{_TORS}
<vdW version="0.4" potential="Lennard-Jones-12-6" combining_rules="Lorentz-Berthelot"
 scale12="{vdw_scale12}" scale13="{vdw_scale13}" scale14="{vdw_scale14}" scale15="{vdw_scale15}"
 cutoff="9.0*angstrom" switch_width="1.0*angstrom"
 periodic_method="cutoff" nonperiodic_method="no-cutoff">{vdw}</vdW>
{charge_block}
<VirtualSites version="0.3" exclusion_policy="{exclusion_policy}">{vs}</VirtualSites>
<Electrostatics version="0.4" scale12="{elec_scale12}" scale13="{elec_scale13}"
 scale14="{elec_scale14}" scale15="{elec_scale15}"
 cutoff="9.0*angstrom" switch_width="0.0*angstrom"
 periodic_potential="Ewald3D-ConductingBoundary" nonperiodic_potential="Coulomb"
 exception_potential="Coulomb"></Electrostatics>
</SMIRNOFF>
"""


def via_interchange(offxml_or_ff, topology, *, charge_from_molecules=None,
                    combine_nonbonded_forces=True):
    """Pathway 1: SMIRNOFF -> Interchange -> OpenMM System.

    Goes through ``create_interchange`` + ``to_openmm_system`` so the
    ``combine_nonbonded_forces`` layout can be controlled (needed to exercise
    interchange #1498).
    """
    ff = offxml_or_ff if isinstance(offxml_or_ff, OFFForceField) else OFFForceField(offxml_or_ff)
    interchange = ff.create_interchange(topology, charge_from_molecules=charge_from_molecules)
    return interchange.to_openmm_system(combine_nonbonded_forces=combine_nonbonded_forces)


def via_openmmforcefields(molecules, offxml, omm_topology, positions, *,
                          nonbonded_method=NoCutoff, preload=True, **create_kwargs):
    """Pathway 2: openmmforcefields SMIRNOFFTemplateGenerator + createSystem.

    Returns ``(system, positions_with_vsites)``.  ``preload=True`` generates and
    loads every residue template before ``createSystem`` (belt-and-suspenders
    against a historical OpenMM constraint-leak bug, fixed in OpenMM PR #5236 /
    8.5.x -- harmless when already fixed).
    """
    if isinstance(molecules, (list, tuple)):
        mols = list(molecules)
    else:
        mols = [molecules]
    generator = SMIRNOFFTemplateGenerator(molecules=mols, forcefield=offxml)
    omm_ff = OMMForceField()
    omm_ff.registerTemplateGenerator(generator.generator)
    if preload:
        # Dedupe isomorphic molecules so we don't try to register the same
        # residue template twice (which raises "already exists").
        loaded = []
        for mol in mols:
            if any(mol.is_isomorphic_with(seen) for seen in loaded):
                continue
            omm_ff.loadFile(io.StringIO(generator.generate_residue_template(mol)))
            loaded.append(mol)
    modeller = Modeller(omm_topology, positions)
    modeller.addExtraParticles(omm_ff)
    system = omm_ff.createSystem(modeller.topology, nonbondedMethod=nonbonded_method,
                                 **create_kwargs)
    return system, modeller.positions
