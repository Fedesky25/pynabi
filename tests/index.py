from pynabi import Vec3D, Atom, DataSet, AbIO, I, O, createAbi, SymmetricKGrid, BZ, KPath, presets, SKO, SpinType

Si = Atom.of("Si")
Zn = Atom.of("Zn")
S = Atom.of("S")

base = DataSet(
    presets.ZincBlende(5.17, Zn, S, "./pseudos/PBE"),
    AbIO().prefix(output="scf/bande")
)

d1 = DataSet(
    SymmetricKGrid(
        BZ.Irreducible, 4,
        Vec3D(0.5, 0.5, 0.5),
        Vec3D(0.5, 0.0, 0.0),
        Vec3D(0.0, 0.5, 0.0),
        Vec3D(0.0, 0.0, 0.5),
    ),
    AbIO().print(O.ElectronDensity)
)
d2 = DataSet(
    KPath(10, 
            Vec3D.zero(), 
            Vec3D(0.5, 0, 0.5),
            Vec3D(0.5, 0, 0.5),
            Vec3D(0.5, 0.25, 0.75),
            Vec3D(0.5, 0.5, 0.5),
            Vec3D.zero(),
            Vec3D(3/8, 3/8, 3/4),
            Vec3D(1/2, 0, 1/2)
    ),
    AbIO().get(I.ElectronDensity, d1)
)

txt = createAbi(base, d1, d2)
with open("./prova.txt", "w") as f:
    f.write(txt)
