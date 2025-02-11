import math
from typing import Sequence

import openmm
import openmm.app
import openmm.unit
from openff.interchange import Interchange
from openff.interchange.components.potentials import Potential
from openff.interchange.models import (
    LibraryChargeTopologyKey,
    PotentialKey,
    SingleAtomChargeTopologyKey,
    TopologyKey,
)
from openff.toolkit import ForceField, Molecule, Quantity, Topology
from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper


def get_charge_sum(
    interchange: Interchange,
    topology_indices: Sequence[int],
) -> Quantity:
    charges = {
        key.atom_indices[0]: value.m
        for key, value in interchange["Electrostatics"].charges.items()
    }

    return Quantity(
        sum([charges[index] for index in topology_indices]),
        "elementary_charge",
    )


def smear_charges(
    interchange: Interchange,
    topology_indices: Sequence[int],
    charge_to_smear: Quantity,
) -> Interchange:
    initial_charge_sum = get_charge_sum(interchange, topology_indices)

    per_atom_difference = charge_to_smear / len(topology_indices)

    interchange["Electrostatics"]._charges_cached = False

    for index in topology_indices:
        topology_key = SingleAtomChargeTopologyKey(this_atom_index=index)
        potential_key = interchange["Electrostatics"].key_map[topology_key]

        interchange["Electrostatics"].potentials[potential_key].parameters[
            "charge"
        ] -= per_atom_difference

    new_charge_sum = get_charge_sum(interchange, topology_indices)

    assert math.isclose((initial_charge_sum - new_charge_sum).m, charge_to_smear.m)

    return interchange


def get_total_charge(system: openmm.System) -> float:
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            return sum(
                [
                    force.getParticleParameters(index)[0]._value
                    for index in range(force.getNumParticles())
                ]
            )


def run(simulation: openmm.app.Simulation):
    """Run an OpenMM simulation for 1 minutes and write a trajectory."""

    dcd_reporter = openmm.app.DCDReporter("trajectory.dcd", 100)
    simulation.reporters.append(dcd_reporter)

    simulation.context.computeVirtualSites()
    simulation.minimizeEnergy()
    simulation.context.setVelocitiesToTemperature(
        simulation.integrator.getTemperature()
    )
    simulation.runForClockTime(1.0 * openmm.unit.minute)


sage_ff14sb = ForceField("openff-2.2.0.offxml", "ff14sb_off_impropers_0.0.4.offxml")

# Add a dummy 0.0 library charge at the _top_ so it's only used as a last resort
sage_ff14sb["LibraryCharges"].add_parameter(
    parameter_kwargs={
        "smirks": "[*:1]",
        "charge1": Quantity(0.0, "elementary_charge"),
        "name": "dummy",
    },
    before=0,  # "[#3+1:1]",
)

substructure_mol = Molecule.from_smiles(
    "C1=CC2=C(C=C1N3C(=O)C[C@H](SC[C@H](N[Fr])C(=O)[Fr])C3=O)C(=O)O[C@]24C5=C(C=C(C=C5)O)OC6=C4C=CC(=C6)O"
)

for atom in substructure_mol.atoms:
    if atom.symbol == "Fr":
        atom.metadata["substructure_atom"] = False
    else:
        atom.metadata["substructure_atom"] = True

print("loading topoogy ...")
topology = Topology.from_pdb(
    "3ip9_dye_solvated.pdb",
    _additional_substructures=[substructure_mol],
)

protein = topology.molecule(0)

for atom in substructure_mol.atoms:
    if atom.symbol == "Fr":
        atom._atomic_number = 1

print("assigning graph charges ...")
protein.assign_partial_charges(
    partial_charge_method="openff-gnn-am1bcc-0.1.0-rc.3.pt",
    toolkit_registry=NAGLToolkitWrapper(),
)

print("making Interchange ...")
interchange = sage_ff14sb.create_interchange(
    topology,
    charge_from_molecules=[],  # =[protein],
)

potential_keys_to_remove = list()
topology_keys_to_remove = list()

nagl_indices = tuple(
    key.atom_indices[0]
    for key, val in interchange["Electrostatics"].key_map.items()
    if val.id == "[*:1]"
)

print("replacing dummy charges with NAGL charges ... ")
for key, charge in interchange["Electrostatics"].charges.items():
    if key.atom_indices[0] in nagl_indices:
        # only modify charges where dummy placeholder of 0.0 was assigned from "[*:1]" parameter
        index = key.atom_indices[0]

        # must make new "single atom" topology key for each atom since the current
        # 1:many representation from the dummy charge is no longer valid
        new_potential_key = PotentialKey(
            id="inserted_graph_charges",
            associated_handler="molecules_with_preset_charges",
            mult=index,
        )
        new_potential = Potential(parameters={"charge": protein.partial_charges[index]})

        interchange["Electrostatics"].key_map[
            SingleAtomChargeTopologyKey(this_atom_index=index)
        ] = new_potential_key
        interchange["Electrostatics"].potentials.update(
            {new_potential_key: new_potential}
        )

        # remove the keys associated with the dummy library charge
        potential_keys_to_remove.append(
            interchange["Electrostatics"].key_map[
                LibraryChargeTopologyKey(this_atom_index=index)
            ]
        )

        topology_keys_to_remove.append(LibraryChargeTopologyKey(this_atom_index=index))

for key_to_remove in topology_keys_to_remove:
    interchange["Electrostatics"].key_map.pop(key_to_remove)

interchange["Electrostatics"]._charges_cached = False

formal_charge_sum_of_nagl_indices = sum(
    [atom.formal_charge for atom in substructure_mol.atoms]
)

charge_to_smear = (
    get_charge_sum(interchange, nagl_indices) - formal_charge_sum_of_nagl_indices
)

interchange = smear_charges(
    interchange=interchange,
    topology_indices=nagl_indices,
    charge_to_smear=charge_to_smear,
)

assert math.isclose(
    get_charge_sum(interchange, nagl_indices).m,
    formal_charge_sum_of_nagl_indices.m,
    abs_tol=1e-10,
    rel_tol=0,
)

print("making OpenMM simulation ...")
simulation = interchange.to_openmm_simulation(
    integrator=openmm.LangevinMiddleIntegrator(
        300 * openmm.unit.kelvin,
        1 / openmm.unit.picosecond,
        0.004 * openmm.unit.picoseconds,
    ),
)

print(f"total system charge is {get_total_charge(simulation.system)}")

print("serializing OpenMM system ...")
with open("system.xml", "w") as f:
    f.write(openmm.XmlSerializer.serialize(simulation.system))

print("running for 1 minute ...")
run(simulation)
