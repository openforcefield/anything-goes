"""Stress test: 5 alkane types x 10 copies = 50 shuffled — tests template reuse."""
import time, copy, random
import openmm, openmm.unit
from openff.toolkit import ForceField as OFFForceField, Molecule, Topology
from openmm.app import ForceField, Modeller, NoCutoff
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openmm import Vec3

EU = openmm.unit.kilocalories_per_mole
FF = "openff_no_water-3.0.0-alpha0.offxml"
off_ff = OFFForceField(FF)


def energy(sys, pos):
    sys = copy.deepcopy(sys)
    ig = openmm.VerletIntegrator(0.001)
    ctx = openmm.Context(sys, ig, openmm.Platform.getPlatformByName("Reference"))
    ctx.setPositions(pos)
    ctx.applyConstraints(ig.getConstraintTolerance())
    e = ctx.getState(getEnergy=True).getPotentialEnergy()
    del ctx, ig
    return e


def build_top(mols, sp=3.0):
    top = mols[0].to_topology().to_openmm()
    pos = mols[0].conformers[0].to_openmm()
    mod = Modeller(top, pos)
    for i, m in enumerate(mols[1:], 1):
        raw = m.conformers[0].to_openmm().value_in_unit(openmm.unit.nanometer)
        mod.add(
            m.to_topology().to_openmm(),
            [Vec3(x + sp * i, y, z) for x, y, z in raw] * openmm.unit.nanometer,
        )
    return mod


random.seed(42)
uniq = []
for smi in ["CCC", "CCCCC", "CC(C)C", "CC(C)CC", "CCCCCCC"]:
    m = Molecule.from_smiles(smi)
    m.generate_conformers(n_conformers=1)
    uniq.append(m)

full = uniq * 10
random.shuffle(full)
mod = build_top(full)

g = SMIRNOFFTemplateGenerator(molecules=uniq, forcefield=FF)
f = ForceField()
f.registerTemplateGenerator(g.generator)
t0 = time.perf_counter()
of = f.createSystem(mod.topology, nonbondedMethod=NoCutoff)
t_of = time.perf_counter() - t0

ct = Topology()
for m in full:
    ct.add_molecule(m)
t0 = time.perf_counter()
ic = off_ff.create_openmm_system(ct)
t_ic = time.perf_counter() - t0

de = abs((energy(ic, mod.positions) - energy(of, mod.positions)).value_in_unit(EU))
c_ok = ic.getNumConstraints() == of.getNumConstraints()

print(
    f"TEST 2 (5types x 10copies): IC={t_ic:.1f}s OF={t_of:.1f}s "
    f"dE={de:.1e} constraints={'PASS' if c_ok else 'FAIL'} "
    f"atoms={of.getNumParticles()}"
)
