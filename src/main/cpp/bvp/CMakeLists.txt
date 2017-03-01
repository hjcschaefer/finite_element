FILE(GLOB bvp_src *.cpp *.hpp)
ADD_LIBRARY(bvp ${bvp_src})
TARGET_LINK_LIBRARIES(bvp felib)
#TARGET_COMPILE_FEATURES(bvp PRIVATE cxx_range_for)

