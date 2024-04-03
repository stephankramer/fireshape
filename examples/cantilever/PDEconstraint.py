from firedrake import *
from fireshape import PdeConstraint


class LinearElasticitySolver(PdeConstraint):
    """Linear elasticity equation. """

    def __init__(self, mesh_m):
        super().__init__()
        self.mesh_m = mesh_m

        # Setup problem, linear finite elements
        self.V = VectorFunctionSpace(mesh_m, "CG", 1)

        # Preallocate soln variables for state eqn
        self.solution = Function(self.V, name="State")
        self.testfunction = TestFunction(self.V)

        # Define surface dead load and Lame' parameters
        self.g = as_vector([0, -1])
        E = 15
        nu = 0.35
        self.mu = E/(2*(1+nu))
        self.lbd = E*nu/((1+nu)*(1-2*nu))

        # Weak form of linear elasticity
        u = self.solution
        v = self.testfunction
        self.F = (inner(self.sigma(u), self.epsilon(v))*dx +
                  -dot(self.g, v)*ds(2))

        # Dirichlet boundary conditions
        u_fix = as_vector([0, 0])
        self.bcs = DirichletBC(self.V, u_fix, [1])

        # PDE-solver parameters
        self.params = {"ksp_type": "preonly",
                       "pe_type": "lu"}

    def epsilon(self, w):
        return 0.5*(grad(w) + grad(w).T)

    def sigma(self, w):
        Id = Identity(self.mesh_m.geometric_dimension())
        return self.lbd*div(w)*Id + 2*self.mu*self.epsilon(w)

    def solve(self):
        super().solve()
        solve(self.F == 0, self.solution, bcs=self.bcs,
              solver_parameters=self.params)
