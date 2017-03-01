#include <felib/HatFunction.hpp>

// ------- LEFT HALF HAT --------------------------------------
LeftHalfHatFunction::LeftHalfHatFunction(size_t idx, const std::vector<double>& grid) {
    _left = grid[idx];
    _right = grid[idx+1];
}

std::pair<double, double> LeftHalfHatFunction::domain() const {
    return std::make_pair(_left, _right);
}

double LeftHalfHatFunction::value_impl(double x) const {
    return (_right - x) / (_right - _left);
}

double LeftHalfHatFunction::deriv_impl(double x) const {
    return -1.0 / (_right - _left);
}

// ------- RIGHT HALF HAT --------------------------------------
RightHalfHatFunction::RightHalfHatFunction(size_t idx, const std::vector<double>& grid) {
    _left = grid[idx-1];
    _right = grid[idx];
}

std::pair<double, double> RightHalfHatFunction::domain() const {
    return std::make_pair(_left, _right);
}

double RightHalfHatFunction::value_impl(double x) const {
    return (x - _left) / (_right - _left);
}

double RightHalfHatFunction::deriv_impl(double x) const {
    return 1.0 / (_right - _left);
}

bool RightHalfHatFunction::within(double x) const {
    auto om = domain();
    return (x >= om.first) && (x <= om.second);
}
// ------- HAT --------------------------------------

HatFunction::HatFunction(size_t idx, const std::vector<double>& grid) {
    _lower = grid[idx-1];
    _mid = grid[idx];
    _upper = grid[idx+1];
}

std::pair<double, double> HatFunction::domain() const {
    return std::make_pair(_lower, _upper);
}

double HatFunction::value_impl(double x) const {
    if (x < _mid) {
        return (x - _lower) / (_mid - _lower);
    }
    return (_upper - x)/(_upper-_mid);
}

double HatFunction::deriv_impl(double x) const {
    if (x < _mid) {
        return 1.0 / (_mid - _lower);
    }
    return - 1.0/(_upper-_mid);
}

template<>
std::vector<std::unique_ptr<AbstractBasisFunction>> getBasis<HatFunction>(const std::vector<double>& grid,
                                                                          bool dirichiletLeft,
                                                                          bool dirichiletRight) {
    std::vector<std::unique_ptr<AbstractBasisFunction>> basis;
    if (!dirichiletLeft) {
        basis.push_back(std::unique_ptr<AbstractBasisFunction>(new LeftHalfHatFunction(0, grid)));
    }
    for(size_t i=1; i < grid.size()-1; i++) {
        basis.push_back(std::unique_ptr<AbstractBasisFunction>(new HatFunction(i, grid)));
    }
    if (!dirichiletRight) {
        basis.push_back(std::unique_ptr<AbstractBasisFunction>(new RightHalfHatFunction(grid.size()-1, grid)));
    }
    return basis;
}
