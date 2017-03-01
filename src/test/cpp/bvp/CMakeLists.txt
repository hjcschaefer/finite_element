FILE(GLOB bvp_test_src *.cpp)
ADD_EXECUTABLE(bvp-test ${bvp_test_src})
ADD_TEST(bvp-test bvp-test)
TARGET_LINK_LIBRARIES(bvp-test bvp gtest gtest_main)
TARGET_COMPILE_FEATURES(bvp-test PUBLIC cxx_range_for)

