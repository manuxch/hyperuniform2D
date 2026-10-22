cmake_minimum_required(VERSION 3.20)
project(hyperuniform VERSION 2.1.0 LANGUAGES CXX)

# ── Estándar ────────────────────────────────────────────────────────────────
set(CMAKE_CXX_STANDARD          20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS        OFF)   # -std=c++20, no -std=gnu++20

# ── Tipo de build por defecto ────────────────────────────────────────────────
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Tipo de compilación" FORCE)
endif()

# ── Fuentes ─────────────────────────────────────────────────────────────────
# Todo salvo main.cpp va a una biblioteca estática que comparten el ejecutable y los tests.
add_library(hyper2d_core STATIC
    src/parametros.cpp
    src/config.cpp
    src/fourier.cpp
    src/annealing.cpp
    src/io.cpp
    src/thread_pool.cpp
    src/sha256.cpp
    src/obs.cpp
    src/continuum.cpp
    src/minimize.cpp
    src/continuum_driver.cpp
    src/continuum_mc.cpp
)
target_include_directories(hyper2d_core PUBLIC include)
target_compile_definitions(hyper2d_core PUBLIC HYPER2D_VERSION="${PROJECT_VERSION}")

find_package(Threads REQUIRED)
target_link_libraries(hyper2d_core PUBLIC Threads::Threads)

add_executable(hyper2D src/main.cpp)
target_link_libraries(hyper2D PRIVATE hyper2d_core)

add_executable(test_hyper2D tests/test_hyper2D.cpp)
target_link_libraries(test_hyper2D PRIVATE hyper2d_core)

# -march=native optimiza para la CPU donde se compila: en un cluster, compilar en un nodo de cómputo
# o usar -DHYPER2D_NATIVE=OFF si el binario se compila en el nodo de login y corre en otros nodos.
option(HYPER2D_NATIVE "Compilar con -march=native en Release" ON)

# ── Opciones de compilación ──────────────────────────────────────────────────
foreach(tgt hyper2d_core hyper2D test_hyper2D)
    target_compile_options(${tgt} PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang>:
            -Wall -Wextra -Wpedantic
            $<$<CONFIG:Release>:-O3>
            $<$<AND:$<CONFIG:Release>,$<BOOL:${HYPER2D_NATIVE}>>:-march=native>
            $<$<CONFIG:Debug>:-g -fsanitize=address,undefined>
        >
        $<$<CXX_COMPILER_ID:MSVC>:
            /W4
            $<$<CONFIG:Release>:/O2>
        >
    )
    # Address/UB sanitizers también necesitan el flag en el linker
    target_link_options(${tgt} PRIVATE
        $<$<AND:$<CXX_COMPILER_ID:GNU,Clang>,$<CONFIG:Debug>>:
            -fsanitize=address,undefined
        >
    )
endforeach()

# ── Tests ────────────────────────────────────────────────────────────────────
enable_testing()
add_test(NAME hyper2D_unit COMMAND test_hyper2D)
find_package(Python3 COMPONENTS Interpreter)
if(Python3_FOUND)
    add_test(NAME hyper2D_obs_end_to_end
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_obs.py
                     $<TARGET_FILE:hyper2D> ${CMAKE_CURRENT_SOURCE_DIR}/scripts/obs_tools.py)
endif()

# ── Información al final de cmake ────────────────────────────────────────────
message(STATUS "hyper2D    : ${PROJECT_VERSION}")
message(STATUS "Build type : ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ standard: C++${CMAKE_CXX_STANDARD}")
message(STATUS "Compiler   : ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
