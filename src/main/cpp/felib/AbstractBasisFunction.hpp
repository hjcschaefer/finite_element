#pragma once 
#ifndef ABSTRACTBASISFUNCTION_HPP
#define ABSTRACTBASISFUNCTION_HPP 

#include <tuple>
#include <memory>
#include <vector>

class AbstractBasisFunction {
public:
    virtual std::pair<double, double> domain() const = 0;
    double value(double x) const {
        if (!within(x)) return 0.0;
        return value_impl(x);
    }
    double deriv(double x) const {
        if (!within(x)) return 0.0;
        return deriv_impl(x);
    }

private:
    virtual double value_impl(double x) const = 0;
    virtual double deriv_impl(double x) const = 0;

    virtual bool within(double x) const {
        auto om = domain();
        return (x >= om.first) && (x < om.second);
    }
};

template<typename BasisFunction>
std::vector<std::unique_ptr<AbstractBasisFunction>> getBasis(const std::vector<double>& grid,
                                                             bool dirichiletLeft,
                                                             bool dirichiletRight);

#endif /* ABSTRACTBASISFUNCTION_HPP */
