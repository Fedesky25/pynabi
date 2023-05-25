import sys, os
sys.path.append(os.getcwd())

from pynabi import createAbi, DataSet, presets, Atom, BZ, SymmetricKGrid, Vec3D, ToleranceOn, EnergyCutoff, StepNumber, SCFProcedure

pseudo_folder = "./pseudos/PBE-SR"
Si = Atom.of("Si")

base = DataSet(
    presets.BCC(5.09, Si, Si, pseudo_folder),
    SymmetricKGrid(BZ.Irreducible, 2, 
                   Vec3D(0.5, 0.5, 0.5),
                   Vec3D(0.5, 0.0, 0.0),
                   Vec3D(0.0, 0.5, 0.0),
                   Vec3D(0.0, 0.0, 0.5),
                   ),
    ToleranceOn.EnergyDifference(1e-6),
    SCFProcedure(0),
    StepNumber(30)
)
d = [DataSet(EnergyCutoff(8.0 + i*0.25)) for i in range(0,17)]

with open("./out.txt", 'w') as f:
    f.write(createAbi(base, *d))