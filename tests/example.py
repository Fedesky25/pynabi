import sys, os
sys.path.append(os.getcwd())

from pynabi import createAbi, DataSet, AbIn, AbOut, Occupation
from pynabi.kspace import CriticalPointsOf, BZ, SymmetricGrid, path, UsualKShifts
from pynabi.calculation import ToleranceOn, EnergyCutoff, MaxSteps, SCFMixing, NonSelfConsistentCalc
from pynabi.crystal import Atom, FluoriteLike
from pynabi.units import EUnit

# create manually an atom -> Atom(<Z>, <pseudo potential name>)
# or using sensible defaults as follows
Zr = Atom.of("Zr")  # Z=40 and pseudos located at "Zr.psp8"
Oxy = Atom.of("O")  # Z=8 and pseudos located at "O.psp8"

# base dataset with common variables
base = DataSet(
    AbOut("./scf/scf"),                             # prefixes for output files
    AbIn().PseudoPotentials("./pseudos/PBE-SR"),    # folder with pseudo potentials
    FluoriteLike(Zr, Oxy, 5.14),                    # creates AtomBasis and Lattice of a crystal like fluorite
    SymmetricGrid(BZ.Irreducible, UsualKShifts.FCC) # easily define kptopt, ngkpt, nshiftk, kpt
        .ofMonkhorstPack(4),
    SCFMixing(density=True).Pulay(10),              # scf cycle with Pulay mixing of the density based on the last 10 iteration
    ToleranceOn.EnergyDifference(1e-6),             # expressively define the tolerance
    MaxSteps(30)                                    # nstep
)

# set the default energy unit in eV (from now on)
EUnit.eV.setAsDefault()

# datasets to see the convergenge as a function of energy
sets = [DataSet(EnergyCutoff.of(8.0 + i*0.25)) for i in range(0,17)]

# final non-self-consistent round to find bands 
bands = DataSet(
    NonSelfConsistentCalc(),
    ToleranceOn.WavefunctionSquaredResidual(1e-12),
    AbIn().ElectronDensity(sets[-1]),               # get the electron density from the last dataset
    Occupation.EqualBandNumber(2.0, 8),             # same number of bands (max 8) for each k point
    path(10, "GXWKGLUWLK", CriticalPointsOf.FCC)    # easily define a path in the k-space   
)

with open("./out.txt", 'w') as f:
    f.write(createAbi(base, *sets, bands))