__all__ = ["Vec3D", "SKO"]

class Vec3D:
    def __init__(self, x: float, y: float, z: float) -> None:
        self.x = x
        self.y = y
        self.z = z
    
    def __str__(self) -> str:
        return "{} {} {}".format(self.x, self.y, self.z)
    
    @staticmethod
    def uniform(v: float = 1.0):
        return Vec3D(v,v,v)
    
    @staticmethod
    def zero():
        return Vec3D(0,0,0)

class Stampable:
    def stamp(self, index: int):
        raise NotImplementedError(f"{type(self).__name__} has not implemented stamp")
    
    def __str__(self):
        return self.stamp(0)


class OneLineStamp(Stampable):
    name = ''
    
    def __init__(self, value, suffix: str = '') -> None:
        super().__init__()
        self.value = value
        self.suffix = suffix

    def stamp(self, index: int):
        t = type(self)
        return f"{t.name}{self.suffix}{index or ''} {self.value}"


class SKO:
    """Single K Occupation"""
    def __init__(self, k: Vec3D, *occupations: float) -> None:
        self.k = k
        self.occ = occupations


def sectionTitle(index: int, name: str):
    if index > 0:
        return f"\n# DS{index} - {name}"
    else:
        return f"\n# {name}"
    
    
def _pos_int(v):
    return type(v) is int and v > 0