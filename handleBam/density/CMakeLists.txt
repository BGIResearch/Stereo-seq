CMAKE_MINIMUM_REQUIRED(VERSION 3.5)

set(SRC_PATH ${PROJECT_SOURCE_DIR})
#set(CXX_INCLUDE_PATH ${SRC_PATH}/cxx/include)
#set(CXX_LIB_PATH ${SRC_PATH}/cxx/lib)


#set(CMAKE_CXX_FLAGS "-std=c++17 -O2 -DNDEBUG -W -Wall -pedantic -fopenmp -lpthread -lrt")
#set(CMAKE_CXX_FLAGS "-std=c++17 -g -ggdb -gstabs+ -pedantic -fopenmp -lpthread -lrt")


link_directories(${INSTALL_PATH}/lib)

set (src
    kde.cpp
    )

set(LIBRARY_OUTPUT_PATH ${INSTALL_PATH}/lib)
add_library(density ${src})

target_link_libraries(density
    fftw3
    )

set_target_properties(density PROPERTIES OUTPUT_NAME "density")
