# Using the `drivers` module

When using or developing force fields, it can be difficult to know if a force field is being applied correctly - or
even to define what that means.  Molecular dynamics simulations include inherent randomness (barostats, thermostats,
etc.) and have many sources of numerical approximation (PME, float rounding, etc.), so it's not tractable to evaluate
trajectories to validate the application of a force field. Instead, running energy evaluations of systems at the initial
configuration is a high-value, low-effort way to assess correctness, at least correctness in how force field parameters
are applied. Different engines ought to report identical or nearly identical values of energies given the same particle
positions and physics parameters (harmonic force constants, partial charges, etc.). When this is not true, decomposing
the potential energy function allows a viable starting point for debugging. For more on this approach to validating
simulation input files, see [Shirts 2017](https://doi.org/10.1007/s10822-016-9977-1).

Interchange has a high-level API for running these energy evaluations. There are several functions for simple checks
and more options for detailed analyses which have proven useful in developing new features and debugging issues along
the way.

For starters, let's create an octanol-hexane mixture using Interchange's Packmol wrapper:

```python
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
# 140 [34.37913464 34.37913464 34.37913464]
```

We can straightforwardly apply Sage to this topology:

```python
from openff.toolkit import ForceField


interchange = ForceField("openff-2.2.1.offxml").create_interchange(topology)
interchange
# Interchange with 7 collections, periodic topology with 3220 atoms.
```

This `Interchange` object represents the chemical topology - our octanol-hexane mixture - complete with force field
parameters applied, positions of each atom specified, and the periodic box vectors defined. That's all we need to get
single-frame energies from MD engines. (Okay, there are some details in between these steps, but Interchange makes some
reasonable assumptions where necessary.) The `drivers` submodule has functions which make this easy. Each function
internally follows these steps:

1. Export the `Interchange` object to simulation files/objects needed by a particular MD engine
2. Call the engine to run a "zero timestep" simulation
3. Get detailed energy evaluations from the engine
4. Report energies to the user in a tidy `EnergyReport` object

```python
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
```

Each dictionary represents a portion of the overall potential energy function. Numerically, the valence terms
agree to 5-6 significant figures. The non-bonded interactions are a little bit more off, mostly because of minor
differences in how vdW tail corrections are applied and inherent errors in PME calculations. We don't have to squint
to compare them, however, as there is a handy `EnergyReport.compare` method for this, which raises an error if any terms differ too much:

```python
gromacs_energies.compare(openmm_energies)
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
```

The error lists the energy discrepancies for each term whose difference falls outside a configurable tolerance. We can dig into these discrepancies a little more by passing the `detailed=True` argument to each `get_*_energies` function:

```python
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
```

This gives us an even more detailed view into the parts of the potential energy function. We had already split out the
vdW and electrostatics interactions from an overall non-bonded energy, but we we can further split out the 1-4 intramolecular
corrections. This is just a different representation of the same underlying physics - it all adds up to the same
potential energy in the end - but gives us more information:

1. 1-4 interactions show little difference between OpenMM and GROMACS
1. There is probably a bug in the detailed representation of vdW energies from GROMACS, as the vdW and vdW 1-4 terms do not add to the total vdW energy displayed above
1. Differences in non-bonded energies are mostly due to long-range electrostatics
1. Purely by coincidence, the energy of 1-4 electrostatics of this system happen to almost completely cancel out the rest of the electrostatic interactions - -0.4 kJ/mol looks much alongside 10 kJ/mol than +/- 800 kJ/mol!

So far, we've only been comparing the energy evaluations of two engines. This makes it a little difficult to assess correctness since it's not obvious which one should be used as a reference. Interchange helpfully provides a wrapper around these and more functions which also runs these evaluations with Amber and LAMMPS (if installed):

```python
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
get_summary_data(interchange)
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
```

This high-level view makes a few earlier points more obvious:

1. It's (relatively) easy to get valence contributions to match more closely
1. Electrostatic energies are quite difficult to get matching closely
1. There's probably something wrong with how electrostatic energies are processed from Amber (but we think the charges are written correctly!)
1. Different implementations of a vdW switching function / tail correction make have real impacts on energies

So far, we haven't made any modifications to the `Interchange` object since creating it. Why not run an energy minimization and look at the result?

```python
interchange.minimize()

get_summary_data(interchange)
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
```

And there we have it! Energies decreased, and they did so somewhat significantly. Each molecule's internal geometry was
determined by conformer generation with RDKit or OpenEye and each molecule's position in the simulation box was
determined by Packmol - each of these methods are good starting points but not quite at an energy minimum, so it's not
surprising that our energies went down so much.

In this post, we've shown how to use Interchange's `driver` submodule to

1. Get quick high-level energy evaluations from one (or many) simulation engines
1. Compare energy evaluations across different engines (and their respective exports from Interchange)
1. Opt in to detailed energy evaluations
1. Find bugs in the Interchange code itself

What other features might be useful to you? Please [get in touch](https://github.com/openforcefield/openff-interchange/issues) and let us know!