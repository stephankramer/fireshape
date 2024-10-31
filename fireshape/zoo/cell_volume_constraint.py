import firedrake as fd
import fireshape as fs


__all__ = ["MoYoCellVolumeConstraint"]


class MoYoCellVolumeConstraint(fs.ShapeObjective):

    def __init__(self, c, *args, minimum_volume=None, verbose=False, **kwargs):
        super().__init__(*args, **kwargs)
        self.T = self.Q.T
        self.minimum_volume = minimum_volume
        self.c = c
        self.ref_volume = fd.assemble(fd.Constant(1.0) * fd.dx(domain=self.Q.mesh_r))
        self.f = fd.File('test.pvd')
        self.ff = fd.Function(fd.FunctionSpace(self.mesh_m, "DG", 0))

    def value_form(self):
        self.ff.interpolate(fd.max_value(self.minimum_volume - fd.CellVolume(self.mesh_m), 0) / self.minimum_volume)
        self.f.write(self.ff)
        penalty = fd.max_value(self.minimum_volume - fd.CellVolume(self.mesh_m), 0) / self.minimum_volume
        penalty = -fd.ln(fd.max(fd.CellVolume(self.mesh_m) - self.minimum_volume, 1e-10))
        return 1/self.c * penalty / self.ref_volume * fd.dx(domain=self.mesh_m)

    def derivative_form(self, tst):
        return -1/self.c / self.ref_volume * (
                fd.CellVolume(self.mesh_m)/(fd.CellVolume(self.mesh_m) - self.minimum_volume) 
                + fd.ln(fd.CellVolume(self.mesh_m) - self.minimum_volume)) * fd.div(tst) * fd.dx
        return 1/self.c * fd.conditional(fd.gt(self.minimum_volume - fd.CellVolume(self.mesh_m), 0),
                                         (self.minimum_volume - 2*fd.CellVolume(self.mesh_m)) * fd.div(tst),
                                         fd.dot(fd.as_vector((0,0)),tst)
                                         ) / self.ref_volume / self.minimum_volume * fd.dx
