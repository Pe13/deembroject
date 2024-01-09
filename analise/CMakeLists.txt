

add_executable(analise analise.cpp)
target_include_directories(analise PRIVATE include)

target_link_libraries(analise PRIVATE ROOT::Core ROOT::Hist ROOT::Gpad)