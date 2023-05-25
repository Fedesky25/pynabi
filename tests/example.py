import sys, os
sys.path.append(os.getcwd())

from pynabi import createAbi, DataSet, presets, Atom, BZ, SymmetricKGrid, KPath, CriticalPointsOf, Vec3D, ToleranceOn, EnergyCutoff, StepNumber, SCFProcedure, AbIO, I, Occupation

# folder with pseudo potentials
pseudo_folder = "./pseudos/PBE-SR"

# create manually an atom 
# or using sensible defaults as follows
Si = Atom.of("Si") # Z=14 and pseudos located at "Si.psp8"

# base dataset with common variables
base = DataSet(
    AbIO().prefix(output="./scf/scf"),          # prefixes for input, output, temporary files
    presets.BCC(5.09, Si, Si, pseudo_folder),   # creates AtomBasis and Lattice of a BCC
    SymmetricKGrid(BZ.Irreducible, 2,           # easily define kptopt, ngkpt, nshiftk, kpt
                   Vec3D(0.5, 0.5, 0.5),
                   Vec3D(0.5, 0.0, 0.0),
                   Vec3D(0.0, 0.5, 0.0),
                   Vec3D(0.0, 0.0, 0.5)),
    ToleranceOn.EnergyDifference(1e-6),         # expressively define the tolerance
    SCFProcedure(0),                            # iscf
    StepNumber(30)                              # nstep
)

# datasets to see the convergenge as a function of energy
sets = [DataSet(EnergyCutoff(8.0 + i*0.25)) for i in range(0,17)]

# final non-self-consistent round to find bands 
bands = DataSet(
    SCFProcedure(-2),
    ToleranceOn.WavefunctionSquaredResidual(1e-12),
    AbIO().get(I.ElectronDensity, sets[-1]),            # get the electron density from the last dataset
    Occupation.Semiconductor(bands=8),                  # semiconductor occupation (occopt=1) with 8 bands
    KPath.through(10, CriticalPointsOf.BCC, "GHNGPH"),  # define a path in the k-space   
)

with open("./out.txt", 'w') as f:
    f.write(createAbi(base, *sets, bands))