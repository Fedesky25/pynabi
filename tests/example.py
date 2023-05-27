import sys, os
sys.path.append(os.getcwd())

from pynabi import createAbi, DataSet, AbIn, AbOut, presets, Atom, ToleranceOn, EnergyCutoff, StepNumber, SCFProcedure, Occupation
from pynabi.kspace import CriticalPointsOf, BZ, SymmetricGrid, path, UsualKShifts

# folder with pseudo potentials
pseudo_folder = "./pseudos/PBE-SR"

# create manually an atom 
# or using sensible defaults as follows
Si = Atom.of("Si") # Z=14 and pseudos located at "Si.psp8"

# base dataset with common variables
base = DataSet(
    AbOut("./scf/scf"),                             # prefixes for output files
    presets.BCC(5.09, Si, Si, pseudo_folder),       # creates AtomBasis and Lattice of a BCC
    SymmetricGrid(BZ.Irreducible, UsualKShifts.BCC) # easily define kptopt, ngkpt, nshiftk, kpt
        .ofMonkhorstPack(4),
    ToleranceOn.EnergyDifference(1e-6),             # expressively define the tolerance
    SCFProcedure(0),                                # iscf
    StepNumber(30)                                  # nstep
)

# datasets to see the convergenge as a function of energy
sets = [DataSet(EnergyCutoff(8.0 + i*0.25)) for i in range(0,17)]

# final non-self-consistent round to find bands 
bands = DataSet(
    SCFProcedure(-2),
    ToleranceOn.WavefunctionSquaredResidual(1e-12),
    AbIn().ElectronDensity(sets[-1]),                   # get the electron density from the last dataset
    Occupation.Semiconductor(bands=8),                  # semiconductor occupation (occopt=1) with 8 bands
    path(10, "GHNGPH", CriticalPointsOf.BCC),           # easily define a path in the k-space   
)

with open("./out.txt", 'w') as f:
    f.write(createAbi(base, *sets, bands))