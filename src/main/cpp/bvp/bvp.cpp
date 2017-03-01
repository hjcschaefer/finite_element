
/*
 * Simple BVP Problem
 *
 * u'' + sin(\pi x) = 0     u(0) = u(1) = 0
 *
 */

#include <bvp/bvp.hpp>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include <felib/AbstractBasisFunction.hpp>
#include <cmath>

#include <iostream>

//using loadFunctionPtr=double(*)(double x, void* p);

double loadFunction(double x, void *p) {
    AbstractBasisFunction* phi = static_cast<AbstractBasisFunction*>(p);
    return std::sin(M_PI*x)*phi->value(x);
}

double matrixElement(double x, void *p) {
    std::pair<AbstractBasisFunction*,AbstractBasisFunction*>* phis = static_cast<std::pair<AbstractBasisFunction*,AbstractBasisFunction*>*>(p);
    return phis->first->deriv(x)*phis->second->deriv(x);
}


std::vector<double> bvp::solve(const std::vector<double>& grid,
                               const std::vector<std::unique_ptr<AbstractBasisFunction>>& basis) {

    
    size_t dof = basis.size();
    // construct the load vector
    size_t quadSteps = 100;
    gsl_integration_workspace* work = gsl_integration_workspace_alloc(quadSteps);
    gsl_vector* load = gsl_vector_alloc(dof);
    for(size_t i=0; i<dof; i++) {
        double sum, err;
        gsl_function f;
        f.function = &loadFunction;
        f.params = static_cast<void*>(basis[i].get());
        auto status = gsl_integration_qag(&f, grid[0], grid[grid.size()-1], 0.001, 0.001, quadSteps, 1, work, &sum, &err);
        gsl_vector_set(load, i, sum);
    }

    // build the stiffness matrix
    // NOTE: we integrate over the full range, should cut it down to just the element
    gsl_matrix* stiff = gsl_matrix_alloc(dof, dof);
    for(size_t i=0; i<dof; i++) {
        for(size_t j=0; j<dof; j++) {
            auto phis = std::make_pair(basis[i].get(), basis[j].get());
            double sum, err;
            gsl_function f;
            f.function = &matrixElement;
            f.params = &(phis);
            auto status = gsl_integration_qag(&f, grid[0], grid[grid.size()-1], 0.001, 0.001, quadSteps, 1, work, &sum, &err);
            gsl_matrix_set(stiff, i, j, sum);
        }
    }

    // solve the system
    gsl_permutation* perm = gsl_permutation_alloc(dof);

    int sig;
    gsl_linalg_LU_decomp(stiff, perm , &sig);

    gsl_vector* res = gsl_vector_alloc(dof);
    gsl_linalg_LU_solve(stiff, perm, load, res);
    // result is in res
    std::vector<double> solution;
    for(size_t i=0; i<dof; i++) {
        solution.push_back(gsl_vector_get(res, i));
    }
    return solution;
}
