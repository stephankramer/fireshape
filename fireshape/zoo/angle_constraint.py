import firedrake as fd
import fireshape as fs


__all__ = ["MoYoAngleConstraint"]


class MoYoAngleConstraint(fs.ShapeObjective):

    def __init__(self, c, bids, *args, max_angle=None, verbose=False, **kwargs):
        super().__init__(*args, **kwargs)
        self.T = self.Q.T
        self.lam = fd.Function(self.T.function_space())
        self.max_angle = max_angle
        self.bids = bids
        self.c = c

        self.nnode = fd.Function(self.V_m, name="node_normals")
        tst = fd.TestFunction(self.V_m)
        n = fd.FacetNormal(self.mesh_m)
        l = fd.FacetArea(self.mesh_m)
        self.nnode_expr = fd.dot(tst, n)/l * fd.ds
        self.verbose = verbose
        if verbose:
            self.anglef = fd.File('angle.pvd')
            P1_m = fd.FunctionSpace(self.mesh_m, "CG", 1)
            self.cosa_f = fd.Function(P1_m, name="cosa")
            self.penalty_f = fd.Function(P1_m, name="penalty")

    def value_form(self):
        self.nnode = fd.assemble(self.nnode_expr, tensor=self.nnode, form_compiler_parameters=self.params)
        cosa = 2 * fd.dot(self.nnode, self.nnode) - 1
        penalty = fd.max_value(fd.cos(self.max_angle) - cosa, 0)
        if self.verbose:
            self.cosa_f.interpolate(cosa)
            self.penalty_f.interpolate(penalty)
            self.anglef.write(self.nnode, self.penalty_f, self.cosa_f)
        return sum(1/self.c * penalty * fd.ds(bid) for bid in self.bids)

    def derivative(self, out):
        super().derivative(out)
        X = fd.SpatialCoordinate(self.mesh_m)
        # we want d/dX [ dot(tst, n)/l*det(dX/dxi) ] dxi_s
        # we get  1/l d/dX [ dot(tst, n)*det(dX/dxi) ] dxi_s
        # need to add d/dX [1/l] dot(tst,n)*det(dX/dxi) dxi_s
        # we know d/dX [1/l det(dX/dxi)] dxi_s = 0 = d/dX [1/l] dx + 1/l d/dX[det(dX/dxi)] dxi_s
        #                                          = d/dX [1/l] dx + 1/l div(arg2) dx
        # so d/dX [1/l] (arg2) = -1/l div(arg2) ?

        tst = fd.TestFunction(self.V_m)
        n = fd.FacetNormal(self.mesh_m)
        l = fd.FacetArea(self.mesh_m)
        arg2 = fd.TrialFunction(self.V_m)
        ddet = fd.div(arg2) - fd.inner(fd.dot(fd.grad(arg2), n), n)
        dndx = fd.assemble(fd.derivative(self.nnode_expr, X)
                           -1/l * fd.dot(tst,n) * ddet * fd.ds, form_compiler_parameters=self.params).petscmat
        dfdn = fd.assemble(fd.derivative(self.value_form(), self.nnode), form_compiler_parameters=self.params)
        dfdx = fd.Cofunction(self.V_m_dual)
        with dfdn.dat.vec_ro as lvec:
            with dfdx.dat.vec as dfvec:
                dndx.multHermitian(lvec, dfvec)
        out.cofun.dat.data[:] += dfdx.dat.data
