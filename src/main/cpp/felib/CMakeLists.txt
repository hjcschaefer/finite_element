FILE(GLOB felib_src *.hpp *.cpp)
ADD_LIBRARY(felib ${felib_src})
TARGET_LINK_LIBRARIES(felib ${GSL_LIBRARIES})
TARGET_COMPILE_FEATURES(felib PUBLIC cxx_range_for)

