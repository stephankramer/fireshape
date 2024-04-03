from firedrake import *
from fireshape import ShapeObjective
from PDEconstraint import LinearElasticitySolver
import numpy as np


class Compliance(ShapeObjective):
    """Compliance functional for linear Elasticity."""

    def __init__(self,
                 pde_solver: LinearElasticitySolver,
                 Q, cb, *args, **kwargs):
        super().__init__(Q, cb, *args, **kwargs)
        self .pde_solver = pde_solver
        Vdet = FunctionSpace(Q.mesh_r, "DG", 0)
        self.detDT = Function(Vdet)

    def value_form(self):
        """Evaluate compliance functional. """

        #self.detDT.interpolate(det(grad(self.Q.T)))
        #if min(self.detDT.vector()) > 0.05:
        if True:
            u = self.pde_solver.solution
            sigma_u = self.pde_solver.sigma(u)
            eps_u = self.pde_solver.epsilon(u)
            return inner(sigma_u, eps_u)*dx
        else:
            return np.nan * dx(self.pde_solver.mesh_m)
