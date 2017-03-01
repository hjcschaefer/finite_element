#pragma once 
#ifndef HATFUNCTION_HPP
#define HATFUNCTION_HPP 

#include <felib/AbstractBasisFunction.hpp>

#include <vector>

class LeftHalfHatFunction : public AbstractBasisFunction {
public:
    // idx -> idx + 1
    LeftHalfHatFunction(size_t idx, const std::vector<double>& grid);
    virtual std::pair<double, double> domain() const override;

private:
    virtual double value_impl(double x) const override;
    virtual double deriv_impl(double x) const override;

    double _left;
    double _right;
};

class RightHalfHatFunction : public AbstractBasisFunction {
public:
    // idx -1 -> idx
    RightHalfHatFunction(size_t idx, const std::vector<double>& grid);
    virtual std::pair<double, double> domain() const override;

private:
    virtual double value_impl(double x) const override;
    virtual double deriv_impl(double x) const override;
    // overrides from AbstractBasisFunction! need to include rhs;
    virtual bool within(double x) const override;

    double _left;
    double _right;
};

class HatFunction : public AbstractBasisFunction {
public:
    // HatFunction i is centered on _grid[i]
    HatFunction(size_t idx, const std::vector<double>& grid);

    virtual std::pair<double, double> domain() const override;
private:
    virtual double value_impl(double x) const override;
    virtual double deriv_impl(double x) const override;

private:
    double _lower;
    double _mid;
    double _upper;
};


#endif /* HATFUNCTION_HPP */
