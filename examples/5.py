import sys, os
sys.path.append(os.getcwd())

from pynabi import createAbi, DataSet, AbOut, AbIn
from pynabi.crystal import Atom, ZincBlendeLike
from pynabi.calculation import EnergyCutoff, ToleranceOn, SCFMixing, MaxSteps
from pynabi.kspace import SymmetricGrid, BZ, UsualKShifts, path, CriticalPointsOf
from pynabi.occupation import OccupationPerBand


Si = Atom.of("Si")

base = DataSet(
    AbOut("./estruct"),
    AbIn().PseudoPotentials("./pseudos/PBE-SR"),
    ZincBlendeLike(Si, Si, 7.31),
    EnergyCutoff(12.0)
)

scf = DataSet(
    MaxSteps(20),
    SCFMixing().Pulay(10),
    ToleranceOn.EnergyDifference(1e-6),
    SymmetricGrid(BZ.Irreducible, UsualKShifts.FCC)
        .ofMonkhorstPack(4),
)

bands = DataSet(
    AbIn().ElectronDensity(scf),
    ToleranceOn.WavefunctionSquaredResidual(1e-12),
    path(10, "GXWLGKX", CriticalPointsOf.FCC),
    OccupationPerBand(1.0, repeat=8)
)

with open("examples/5.txt", 'w') as file:
    file.write(createAbi(base, scf, bands))