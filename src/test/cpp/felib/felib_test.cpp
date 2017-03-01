#include <gtest/gtest.h>

#include <felib/HatFunction.hpp>

#include <iostream>

TEST(HatFunction, left) {
    std::vector<double> grid {0.0, 0.25, 0.5, 0.75, 1.0};
    LeftHalfHatFunction hat(0, grid);
    ASSERT_NEAR(1.0, hat.value(0.0), 1e-10);
    ASSERT_NEAR(0.5, hat.value(0.25/2), 1e-10);
    ASSERT_NEAR(0.0, hat.value(0.25), 1e-10);
    ASSERT_NEAR(0.0, hat.value(0.35), 1e-10);
}

TEST(HatFunction, right) {
    std::vector<double> grid {0.0, 0.25, 0.5, 0.75, 1.0};
    RightHalfHatFunction hat(grid.size()-1, grid);
    ASSERT_NEAR(1.0, hat.value(1.0), 1e-10);
    ASSERT_NEAR(0.5, hat.value(0.75+0.25/2), 1e-10);
    ASSERT_NEAR(0.0, hat.value(0.75), 1e-10);
    ASSERT_NEAR(0.0, hat.value(0.65), 1e-10);
}

TEST(HatFunction, central) {
    std::vector<double> grid {0.0, 0.25, 0.5, 0.75, 1.0};
    HatFunction hat(1, grid);
    ASSERT_NEAR(1.0, hat.value(0.25), 1e-10);
    ASSERT_NEAR(0.5, hat.value(0.25+0.25/2), 1e-10);
    ASSERT_NEAR(0.5, hat.value(0.25-0.25/2), 1e-10);
    ASSERT_NEAR(0.0, hat.value(0.5), 1e-10);
    ASSERT_NEAR(0.0, hat.value(0.0), 1e-10);
    ASSERT_NEAR(0.0, hat.value(0.65), 1e-10);
}

TEST(BasisFunctionSet, noDirichilet) {
    std::vector<double> grid {0.0, 0.25, 0.5, 0.75, 1.0};
    auto basis = getBasis<HatFunction>(grid, false, false);
    ASSERT_EQ(grid.size(), basis.size());

    double x = 0.0;
    while (x <= 1.0) {
        double sum = 0.0;
        for(auto& p : basis) {
            sum += p->value(x);
        }
        ASSERT_NEAR(1.0, sum, 1e-10);
        x += 0.1;
    }
}

TEST(BasisFunctionSet, dirichilet) {
    std::vector<double> grid {0.0, 0.25, 0.5, 0.75, 1.0};
    auto basis = getBasis<HatFunction>(grid, true, true);
    ASSERT_EQ(grid.size()-2, basis.size());
    double x = 0.25;
    while (x <= 0.75) {
        double sum = 0.0;
        for(auto& p : basis) {
            sum += p->value(x);
        }
        ASSERT_NEAR(1.0, sum, 1e-10);
        x += 0.1;
    }
}


