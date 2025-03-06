import firedrake as fd
import fireshape as fs


__all__ = ["InversionPenalty"]


ip_pvd = fd.VTKFile('InversionPenalty.pvd')
class InversionPenalty(fs.DeformationObjective):
    def __init__(self, min_element_volume, Q, **kwargs):
        super().__init__(Q, **kwargs)
        self.min_element_volume = min_element_volume
        self.P0 = fd.FunctionSpace(Q.mesh_r, "DG", 0)
        self.detJ = fd.Function(self.P0, name='element orientation')
        self.detJ.interpolate(abs(fd.det(fd.Jacobian(Q.mesh_r))))
        self._value_expr = 0.5 * fd.det(fd.grad(self.Q.T)) * self.detJ
        self._value_func = fd.Function(self.P0, name='value')
        self._value_form = fd.ln(0.5 * fd.det(fd.grad(self.Q.T)) * self.detJ - min_element_volume) * fd.dx

    def value_form(self):
        self._value_func.interpolate(self._value_expr)
        ip_pvd.write(self._value_func)
        return self._value_form

    def value(self, x, tol):
        return super().value(x, tol)

    def derivative_form(self, test):
        T = self.Q.T
        return fd.derivative(self.value_form(), T, test)
