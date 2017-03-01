FILE(GLOB felib_test_src *.cpp)
ADD_EXECUTABLE(felib-test ${felib_test_src})
ADD_TEST(felib-test felib-test)
TARGET_LINK_LIBRARIES(felib-test felib gtest gtest_main)
TARGET_COMPILE_FEATURES(felib-test PUBLIC cxx_range_for)

