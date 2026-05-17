cmake_minimum_required(VERSION 3.20)
project(hyperuniform LANGUAGES CXX)

# ── Estándar ────────────────────────────────────────────────────────────────
set(CMAKE_CXX_STANDARD          20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS        OFF)   # -std=c++20, no -std=gnu++20

# ── Tipo de build por defecto ────────────────────────────────────────────────
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Tipo de compilación" FORCE)
endif()

# ── Fuentes ─────────────────────────────────────────────────────────────────
set(SOURCES
    src/main.cpp
    src/parametros.cpp
    src/config.cpp
    src/fourier.cpp
    src/annealing.cpp
    src/io.cpp
    src/thread_pool.cpp
)

add_executable(hyper2D ${SOURCES})

target_include_directories(hyper2D PRIVATE include)

find_package(Threads REQUIRED)
target_link_libraries(hyper2D PRIVATE Threads::Threads)

# ── Opciones de compilación ──────────────────────────────────────────────────
target_compile_options(hyper2D PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:
        -Wall -Wextra -Wpedantic
        $<$<CONFIG:Release>:-O3 -march=native>
        $<$<CONFIG:Debug>:-g -fsanitize=address,undefined>
    >
    $<$<CXX_COMPILER_ID:MSVC>:
        /W4
        $<$<CONFIG:Release>:/O2>
    >
)

# Address/UB sanitizers también necesitan el flag en el linker
target_link_options(hyper2D PRIVATE
    $<$<AND:$<CXX_COMPILER_ID:GNU,Clang>,$<CONFIG:Debug>>:
        -fsanitize=address,undefined
    >
)

# ── Información al final de cmake ────────────────────────────────────────────
message(STATUS "Build type : ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ standard: C++${CMAKE_CXX_STANDARD}")
message(STATUS "Compiler   : ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
