cmake_minimum_required(VERSION 3.5)

project(psKC LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(psKC
    main.cpp

    csr.cpp
    file_name.cpp
    libsvm.cpp

# use by the IDE
    csr.h
    file_name.h
    libsvm.h
    parse_param.h
)

install(TARGETS psKC
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
