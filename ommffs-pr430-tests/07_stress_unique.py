"""Stress test: 20 unique functionalized alkanes — tests many distinct templates."""
import time, copy, itertools, random
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
subs = ["", "(F)", "(N)", "(Cl)"]
cands = []
for bl in [3, 4]:
    for pat in itertools.product(range(len(subs)), repeat=bl):
        cands.append("".join("C" + subs[s] for s in pat))
random.shuffle(cands)

u3 = []
seen = set()
for smi in cands:
    if len(u3) >= 20:
        break
    try:
        m = Molecule.from_smiles(smi, allow_undefined_stereo=True)
        c = m.to_smiles()
        if c in seen:
            continue
        m.generate_conformers(n_conformers=1)
        off_ff.create_openmm_system(m.to_topology())
        seen.add(c)
        u3.append(m)
    except Exception:
        pass

print(f"Validated {len(u3)} unique molecules")
mod = build_top(u3)

g = SMIRNOFFTemplateGenerator(molecules=u3, forcefield=FF)
f = ForceField()
f.registerTemplateGenerator(g.generator)
t0 = time.perf_counter()
of = f.createSystem(mod.topology, nonbondedMethod=NoCutoff)
t_of = time.perf_counter() - t0

ct = Topology()
for m in u3:
    ct.add_molecule(m)
t0 = time.perf_counter()
ic = off_ff.create_openmm_system(ct)
t_ic = time.perf_counter() - t0

de = abs((energy(ic, mod.positions) - energy(of, mod.positions)).value_in_unit(EU))
c_ok = ic.getNumConstraints() == of.getNumConstraints()

print(
    f"TEST 3 ({len(u3)} unique): IC={t_ic:.1f}s OF={t_of:.1f}s "
    f"dE={de:.1e} constraints={'PASS' if c_ok else 'FAIL'} "
    f"atoms={of.getNumParticles()}"
)
