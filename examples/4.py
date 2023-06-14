import sys, os
sys.path.append(os.getcwd())

from pynabi import createAbi, DataSet, AbOut, AbIn
from pynabi.units import Ang, Ha
from pynabi.crystal import Atom, NiAsLike
from pynabi.relaxation import CellOptimization, StructuralOptimization, MaxLatticeDilatation
from pynabi.calculation import EnergyCutoff, ToleranceOn, SCFMixing, MaxSteps
from pynabi.kspace import SymmetricGrid, BZ, UsualKShifts

# Achavalite (FeSe)

Fe = Atom.of("Fe")
Se = Atom.of("Se")

base = DataSet(
    AbOut("./out/FeSe"),
    AbIn().PseudoPotentials("./pseudos/PBE-SR"),
    NiAsLike(Fe, Se, 3.636, 5.946, Ang),
    SymmetricGrid(BZ.Irreducible, UsualKShifts.HEX)
        .ofMonkhorstPack(6),
    EnergyCutoff(16*Ha)
)

relax = DataSet(
    StructuralOptimization(maxSteps=40).BFGS(),
    CellOptimization(0.5*Ha).OnVolumeOnly(),
    MaxLatticeDilatation(5)
)

scf = DataSet(
    SCFMixing().Pulay(),
    MaxSteps(10),
    ToleranceOn.EnergyDifference(1e-6),
    AbIn().WavefunctionsK(relax),
    AbOut().ElectronDensity()
)

with open("examples/4.txt", 'w') as file:
    file.write(createAbi(base, relax, scf))