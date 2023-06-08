import sys, os
sys.path.append(os.getcwd())

from pynabi import createAbi, DataSet, AbOut, AbIn
from pynabi.crystal import Atom, AtomBasis, Lattice
from pynabi.calculation import EnergyCutoff, ToleranceOn, MaxSteps, NonSelfConsistentCalc
from pynabi.kspace import SymmetricGrid, BZ, UsualKShifts
from pynabi.occupation import Metal, Smearing


Si = Atom.of("Si")
Al = Atom.of("Al")

base = DataSet(
    AbOut("./estruct"),
    AbIn().PseudoPotentials("./pseudos/PBE-SR"),
    AtomBasis.ofOne(Atom.of("Al")),
    Lattice.FCC(7.626),
    EnergyCutoff(14.0),
    Metal(Smearing.Marzari5634, 0.04)
)

scf = DataSet(
    MaxSteps(20),
    AbOut().ElectronDensity().ElectronLocalizedFunction(),
    ToleranceOn.EnergyDifference(1e-6),
    SymmetricGrid(BZ.Irreducible, UsualKShifts.FCC)
        .ofMonkhorstPack(6)
)

bands = DataSet(
    AbIn().ElectronDensity(scf),
    AbOut().DensityOfStates(2).FermiSurface(),
    SymmetricGrid(BZ.Irreducible, UsualKShifts.Unshifted)
        .ofMonkhorstPack(30),
    NonSelfConsistentCalc(-3),
    ToleranceOn.WavefunctionSquaredResidual(1e-12),
    Metal.setBands(6),
)


with open("examples/6.txt", 'w') as file:
    file.write(createAbi(base, scf, bands))