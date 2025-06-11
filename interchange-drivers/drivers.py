from rich.pretty import pprint

from openff.toolkit import Molecule, Quantity
from openff.interchange.components._packmol import pack_box, UNIT_CUBE


octanol = Molecule.from_smiles(8 * "C" + "O")
hexane = Molecule.from_smiles(6 * "C")

for molecule in [octanol, hexane]:
    molecule.generate_conformers(n_conformers=1)

topology = pack_box(
    molecules=[octanol, hexane],
    number_of_copies=[60, 80],
    box_shape=UNIT_CUBE,
    target_density=Quantity("800 kg/m**3")
)

print(topology.n_molecules, topology.box_vectors.m_as("angstrom").diagonal())

from openff.toolkit import ForceField


interchange = ForceField("openff-2.2.1.offxml").create_interchange(topology)
interchange
# Interchange with 7 collections, periodic topology with 3220 atoms.

from openff.interchange.drivers.gromacs import get_gromacs_energies
from openff.interchange.drivers.openmm import get_openmm_energies


gromacs_energies = get_gromacs_energies(interchange)
pprint(gromacs_energies)
# EnergyReport(
# │   energies={
# │   │   'Bond': <Quantity(672.310425, 'kilojoule / mole')>,
# │   │   'Angle': <Quantity(3982.40161, 'kilojoule / mole')>,
# │   │   'Torsion': <Quantity(4809.75391, 'kilojoule / mole')>,
# │   │   'vdW': <Quantity(-269.391632, 'kilojoule / mole')>,
# │   │   'Electrostatics': <Quantity(9.76124573, 'kilojoule / mole')>
# │   }
# )

openmm_energies = get_openmm_energies(interchange, combine_nonbonded_forces=False)
pprint(openmm_energies)
# EnergyReport(
# │   energies={
# │   │   'Bond': <Quantity(672.308639, 'kilojoule / mole')>,
# │   │   'Angle': <Quantity(3982.40428, 'kilojoule / mole')>,
# │   │   'Torsion': <Quantity(4809.75712, 'kilojoule / mole')>,
# │   │   'vdW': <Quantity(-269.551771, 'kilojoule / mole')>,
# │   │   'Electrostatics': <Quantity(10.1593947, 'kilojoule / mole')>
# │   }
# )

try:
    gromacs_energies.compare(openmm_energies)
except Exception as error:
    import ipdb; ipdb.set_trace()
    print(error)

# ---------------------------------------------------------------------------
# EnergyError                               Traceback (most recent call last)
# Cell In[5], line 1
# ----> 1 gromacs_energies.compare(openmm_energies)
# 
# File ~/software/openff-interchange/openff/interchange/drivers/report.py:146, in EnergyReport.compare(self, other, tolerances)
#     143         errors[key] = diff
#     145 if errors:
# --> 146     raise EnergyError(errors)
# 
# EnergyError: {'Bond': <Quantity(0.0017857753, 'kilojoule / mole')>, 'Angle': <Quantity(-0.00267123177, 'kilojoule / mole')>, 'Torsion': <Quantity(-0.00321735958, 'kilojoule / mole')>, 'vdW': <Quantity(0.160138541, 'kilojoule / mole')>, 'Electrostatics': <Quantity(-0.398271038, 'kilojoule / mole')>}

pprint(get_gromacs_energies(interchange, detailed=True))
# EnergyReport(
# │   energies={
# │   │   'Bond': <Quantity(672.310425, 'kilojoule / mole')>,
# │   │   'Angle': <Quantity(3982.40161, 'kilojoule / mole')>,
# │   │   'Torsion': <Quantity(4809.75391, 'kilojoule / mole')>,
# │   │   'vdW': <Quantity(-1868.99158, 'kilojoule / mole')>,
# │   │   'vdW 1-4': <Quantity(2059.97388, 'kilojoule / mole')>,
# │   │   'Electrostatics': <Quantity(818.116837, 'kilojoule / mole')>,
# │   │   'Electrostatics 1-4': <Quantity(-808.355591, 'kilojoule / mole')>
# │   }
# )

pprint(get_openmm_energies(interchange, combine_nonbonded_forces=False, detailed=True))
# EnergyReport(
# │   energies={
# │   │   'Bond': <Quantity(672.308639, 'kilojoule / mole')>,
# │   │   'Angle': <Quantity(3982.40428, 'kilojoule / mole')>,
# │   │   'Torsion': <Quantity(4809.75712, 'kilojoule / mole')>,
# │   │   'vdW': <Quantity(-2329.52674, 'kilojoule / mole')>,
# │   │   'vdW 1-4': <Quantity(2059.97497, 'kilojoule / mole')>,
# │   │   'Electrostatics': <Quantity(818.515594, 'kilojoule / mole')>,
# │   │   'Electrostatics 1-4': <Quantity(-808.356199, 'kilojoule / mole')>
# │   }
# )

from openff.interchange.drivers.all import get_summary_data, get_all_energies


pprint(get_all_energies(interchange))
# /Users/mattthompson/software/openff-interchange/openff/interchange/components/mdconfig.py:401: SwitchingFunctionNotImplementedWarning: A switching distance 8.0 angstrom was specified by the force field, but Amber does not implement a switching function. Using a hard cut-off instead. Non-bonded interactions will be affected.
#   warnings.warn(
# /Users/mattthompson/software/openff-interchange/openff/interchange/components/mdconfig.py:277: SwitchingFunctionNotImplementedWarning: A switching distance 8.0 angstrom was specified by the force field, but LAMMPS may not implement a switching function as specified by SMIRNOFF. Using a hard cut-off instead. Non-bonded interactions will be affected.
#   warnings.warn(
#
# {
# │   'OpenMM': EnergyReport(
# │   │   energies={
# │   │   │   'Bond': <Quantity(672.308639, 'kilojoule / mole')>,
# │   │   │   'Angle': <Quantity(3982.40428, 'kilojoule / mole')>,
# │   │   │   'Torsion': <Quantity(4809.75712, 'kilojoule / mole')>,
# │   │   │   'vdW': <Quantity(-269.551771, 'kilojoule / mole')>,
# │   │   │   'Electrostatics': <Quantity(10.1593947, 'kilojoule / mole')>
# │   │   }
# │   ),
# │   'Amber': EnergyReport(
# │   │   energies={
# │   │   │   'Bond': <Quantity(672.30855, 'kilojoule / mole')>,
# │   │   │   'Angle': <Quantity(3982.40442, 'kilojoule / mole')>,
# │   │   │   'Torsion': <Quantity(4809.75695, 'kilojoule / mole')>,
# │   │   │   'vdW': <Quantity(-273.205577, 'kilojoule / mole')>,
# │   │   │   'Electrostatics': <Quantity(39226.4933, 'kilojoule / mole')>
# │   │   }
# │   ),
# │   'GROMACS': EnergyReport(
# │   │   energies={
# │   │   │   'Bond': <Quantity(672.310425, 'kilojoule / mole')>,
# │   │   │   'Angle': <Quantity(3982.40161, 'kilojoule / mole')>,
# │   │   │   'Torsion': <Quantity(4809.75391, 'kilojoule / mole')>,
# │   │   │   'vdW': <Quantity(-269.391632, 'kilojoule / mole')>,
# │   │   │   'Electrostatics': <Quantity(9.76112366, 'kilojoule / mole')>
# │   │   }
# │   ),
# │   'LAMMPS': EnergyReport(
# │   │   energies={
# │   │   │   'Bond': <Quantity(672.309552, 'kilojoule / mole')>,
# │   │   │   'Angle': <Quantity(3982.40422, 'kilojoule / mole')>,
# │   │   │   'Torsion': <Quantity(4809.75719, 'kilojoule / mole')>,
# │   │   │   'vdW': <Quantity(-659.786776, 'kilojoule / mole')>,
# │   │   │   'Electrostatics': <Quantity(8.81260509, 'kilojoule / mole')>
# │   │   }
# │   )
# }

# same thing, just as a table
pprint(get_summary_data(interchange))
# /Users/mattthompson/software/openff-interchange/openff/interchange/components/mdconfig.py:401: SwitchingFunctionNotImplementedWarning: A switching distance 8.0 angstrom was specified by the force field, but Amber does not implement a switching function. Using a hard cut-off instead. Non-bonded interactions will be affected.
#   warnings.warn(
# /Users/mattthompson/software/openff-interchange/openff/interchange/components/mdconfig.py:277: SwitchingFunctionNotImplementedWarning: A switching distance 8.0 angstrom was specified by the force field, but LAMMPS may not implement a switching function as specified by SMIRNOFF. Using a hard cut-off instead. Non-bonded interactions will be affected.
#   warnings.warn(
#
#                Bond        Angle      Torsion         vdW  Electrostatics
# OpenMM   672.308639  3982.404283  4809.757124 -269.551771       10.159395
# Amber    672.308550  3982.404420  4809.756948 -273.205577    39226.493270
# GROMACS  672.310425  3982.401611  4809.753906 -269.391632        9.761124
# LAMMPS   672.309552  3982.404216  4809.757186 -659.786776        8.812605

interchange.minimize()

print(get_summary_data(interchange))
# /Users/mattthompson/software/openff-interchange/openff/interchange/components/mdconfig.py:401: SwitchingFunctionNotImplementedWarning: A switching distance 8.0 angstrom was specified by the force field, but Amber does not implement a switching function. Using a hard cut-off instead. Non-bonded interactions will be affected.
#   warnings.warn(
# /Users/mattthompson/software/openff-interchange/openff/interchange/components/mdconfig.py:277: SwitchingFunctionNotImplementedWarning: A switching distance 8.0 angstrom was specified by the force field, but LAMMPS may not implement a switching function as specified by SMIRNOFF. Using a hard cut-off instead. Non-bonded interactions will be affected.
#   warnings.warn(
#
#                Bond       Angle      Torsion          vdW  Electrostatics
# OpenMM   109.458233  334.696692  4287.825448 -4583.424516      -47.190688
# Amber    109.458042  334.696570  4287.825542 -4586.797027    40220.794510
# GROMACS  109.458092  334.696960  4287.824707 -4583.264130      -47.506271
# LAMMPS   109.458097  334.696915  4287.825411 -4973.377790      -48.666073

