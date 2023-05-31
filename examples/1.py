import sys, os
sys.path.append(os.getcwd())

from pynabi import createAbi, DataSet, AbOut, AbIn
from pynabi.crystal import ZincBlendeLike, Atom
from pynabi.calculation import EnergyCutoff, ToleranceOn, MaxSteps
from pynabi.kspace import SymmetricGrid, BZ, UsualKShifts

Si = Atom.of("Si")

d = DataSet(
    AbOut("./scf/scf"),
    AbIn().PseudoPotentials("./pseudos/PBE-SR"),
    ZincBlendeLike(Si, Si, 7.31),
    SymmetricGrid(BZ.Irreducible,UsualKShifts.FCC)
        .ofMonkhorstPack(2),
    EnergyCutoff.of(12.0),
    ToleranceOn.EnergyDifference(1e-6),
    MaxSteps(10)
)

with open("examples/1.txt", 'w') as file:
    file.write(createAbi(d))