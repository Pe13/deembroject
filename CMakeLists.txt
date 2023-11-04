cmake_minimum_required(VERSION 3.22)
project(deembroject)

# abilita il supporto per i test, tra cui l'opzione BUILD_TESTING usata sotto
include(CTest)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# richiedi l'uso di C++17, senza estensioni non-standard offerte dal compilatore usato
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED OFF)
set(CMAKE_CXX_EXTENSIONS OFF)

# abilita warning
string(APPEND CMAKE_CXX_FLAGS " -Wall -Wextra")

# abilita l'address sanitizer e l'undefined-behaviour sanitizer in debug mode
#string(APPEND CMAKE_CXX_FLAGS_DEBUG " -fsanitize=address,undefined -fno-omit-frame-pointer")
#string(APPEND CMAKE_EXE_LINKER_FLAGS_DEBUG " -fsanitize=address,undefined -fno-omit-frame-pointer")

find_package(ROOT 6.26 COMPONENTS CONFIG REQUIRED)
#find_package(SFML 2.5 COMPONENTS graphics REQUIRED)
#find_package(TGUI 1.0 REQUIRED)
find_package(Boost 1.35 COMPONENTS container REQUIRED)

add_executable(deembroject
        main.cpp
        ParticleType.cpp
        Particle.cpp)

target_include_directories(deembroject PRIVATE include)
target_link_libraries(deembroject PRIVATE ROOT::Core ROOT::Hist ROOT::Physics ROOT::Gpad)

#