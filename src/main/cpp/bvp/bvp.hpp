#pragma once 
#ifndef BVP_HPP
#define BVP_HPP 

#include <felib/AbstractBasisFunction.hpp>

namespace bvp {
    std::vector<double> solve(const std::vector<double>& grid,
                              const std::vector<std::unique_ptr<AbstractBasisFunction>>& basis);
}


#endif /* BVP_HPP */
