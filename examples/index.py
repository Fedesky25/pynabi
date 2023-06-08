import sys, os
sys.path.append(os.getcwd())

from pynabi import createAbi, DataSet, AbIn, AbOut
from pynabi.kspace import CriticalPointsOf, BZ, SymmetricGrid, path, UsualKShifts
from pynabi.calculation import ToleranceOn, EnergyCutoff, MaxSteps, SCFMixing, NonSelfConsistentCalc
from pynabi.crystal import Atom, FluoriteLike
from pynabi.occupation import OccupationPerBand
from pynabi.units import eV, nm

# create manually an atom -> Atom(<Z>, <pseudo potential name>)
# or using sensible defaults as follows
Zr = Atom.of("Zr")  # Z=40 and pseudos located at "Zr.psp8"
Oxy = Atom.of("O")  # Z=8 and pseudos located at "O.psp8"

# base dataset with common variables
base = DataSet(
    AbOut("./scf/scf"),                             # prefix for output files
    AbIn().PseudoPotentials("./pseudos/PBE-SR"),    # folder with pseudo potentials

    FluoriteLike(Zr, Oxy, 0.5135*nm),               # creates AtomBasis and Lattice of a crystal like fluorite
                                                    # with lattice constant 0.5135nm

    SymmetricGrid(BZ.Irreducible, UsualKShifts.FCC) # easily define kptopt, ngkpt, nshiftk, kpt
        .ofMonkhorstPack(4),

    SCFMixing(density=True).Pulay(10),              # scf cycle with Pulay mixing of the density 
                                                    # based on the last 10 iteration

    ToleranceOn.EnergyDifference(1e-6),             # expressively define the tolerance
    MaxSteps(30)                                    # nstep
)

# set the default energy unit in eV (from now on)
eV.setAsReference()

# datasets to see the convergenge as a function of energy
sets = [DataSet(EnergyCutoff(8.0 + i*0.25)) for i in range(0,17)]

# final non-self-consistent round to find bands 
bands = DataSet(
    NonSelfConsistentCalc(),
    ToleranceOn.WavefunctionSquaredResidual(1e-12),
    AbIn().ElectronDensity(sets[-1]),               # get the electron density from the last dataset
    OccupationPerBand(2.0, repeat=8),               # same number of bands (max 8) for each k point
    path(10, "GXWKGLUWLK", CriticalPointsOf.FCC)    # easily define a path in the k-space   
)

with open("./out.txt", 'w') as f:
    f.write(createAbi(base, *sets, bands))