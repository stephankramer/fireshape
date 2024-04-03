from firedrake import *
from firedrake import ConvergenceError

from fireshape import *

import fireshape.zoo as fsz

import ROL

from PDEconstraint import LinearElasticitySolver
from objective import Compliance

# setup problem

mesh = Mesh("cantilever.msh")

def shape_optimize(mesh):
    Q = FeControlSpace(mesh)

    inner = ElasticityInnerProduct(Q, fixed_bids=[1, 2])
    q = ControlVector(Q, inner)

    # setup PDE constraint
    e = LinearElasticitySolver(Q.mesh_m)

    # save state variable evolution in file u.pud
    out = File("solution/u.pvd")




    # create PDE-constrained objective functional
    J_ = None
    def cb():
        out.write(e.solution)
        print("Compliance:", assemble(J_.value_form()))
    J_ = Compliance(e, Q, cb=cb)
    J = ReducedObjective(J_, e)

    # add regularization to improve mesh quality
    Jq = fsz.MoYoSpectralConstraint(100, Constant(0.5), Q)
    Jbbox = fsz.MoYoBoxConstraint(100, (3,), Q, lower_bound=(0,-4), upper_bound=(20,4))
    J = J + Jbbox + Jq

    # Set up volume constraint
    vol = fsz.VolumeFunctional(Q)
    vol0 = vol.value(q, None)

    C = EqualityConstraint([vol], target_value=[vol0])
    M = ROL.StdVector(1)

    # ROL parameters
    params_dict = {
        'General': {'Secant': {'Type': 'Limited-Memory BFGS'}},
        'Step': {'Type': 'Augmented Lagrangian',
                 'Augmented Lagrangian':
                         {'Subproblem Step Type': 'Trust Region',
                          'Subproblem Iteration Limit': 10,
                          'Penalty Parameter Growth Factor': 2}},
        'Status Test': {'Gradient Tolerance': 1e-2,
                        'Step Tolerance': 1e-3,
                        'Constraint Tolerance': 1e-2,
                        'Iteration Limit': 15}}
    params = ROL.ParameterList(params_dict, "Parameters")
    problem = ROL.OptimizationProblem(J, q, econ=C, emul=M)
    solver = ROL.OptimizationSolver(problem, params)
    try:
        solver.solve()
    except ConvergenceError:
        print("ConvergenceError - but let's see if adapt can fix this!")
        pass
    return e.solution

from animate import RiemannianMetric, adapt

fpvd = fd.File('sol.pvd', adaptive=True)
sol = Function(FunctionSpace(mesh, "CG", 1), name='State')
fpvd.write(sol)
for i in range(20):

    sol = shape_optimize(mesh)
    fpvd.write(sol)

    mesh = sol.function_space().mesh()
    TV = fd.TensorFunctionSpace(mesh, "CG", 1)
    metric = RiemannianMetric(TV)
    #metric.compute_hessian(sol[0])
    #metric.enforce_spd()
    #metric2 = RiemannianMetric(TV)
    #metric2.compute_hessian(sol[1])
    #metric2.enforce_spd()
    #metric.intersect(metric2)
    metric.set_parameters({
        'dm_plex_metric_h_min': 0.01,
        'dm_plex_metric_h_max': 0.5,
        'dm_plex_metric_a_max': 5,
        'dm_plex_metric_target_complexity': 100
        })
    metric.normalise()
    mesh = adapt(mesh, metric)
