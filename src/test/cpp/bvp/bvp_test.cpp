#include <gtest/gtest.h>

#include <bvp/bvp.hpp>
#include <felib/HatFunction.hpp>
#include <cmath>

#include <fstream>

TEST(bvp, simple) {
    std::vector<double> grid {0.0, 0.25, 0.5, 0.75, 1.0};
    auto basis = getBasis<HatFunction>(grid, true, true);
    auto res = bvp::solve(grid, basis);

    std::ofstream fout("bla.dat");
    double x = 0.0;
    while (x <= 1.0) {
        double sum = 0.0;
        for(size_t i=0; i<res.size(); i++) {
            sum += res[i]*basis[i]->value(x);
        }
        double ref = sin(M_PI*x)/(M_PI*M_PI);
        std::cerr << x << "\t" << ref << "\t" << sum << std::endl;
        fout << x << "\t" << ref << "\t" << sum << "\t" << sum-ref << std::endl;
        ASSERT_NEAR(ref, sum, 0.007);
        x += 0.05;
    }
    fout.close();

}
