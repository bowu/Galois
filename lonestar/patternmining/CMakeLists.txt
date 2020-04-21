set(sources 
  directed_graphlet.cpp
  execution_plan.cpp
  graphlet.cpp
  hyper_plan.cpp
  multi_plan.cpp
  multip_red_dg.cpp
  multirestplans.cpp
  PatternMining.cpp
  rest_plan.cpp
  rest_set.cpp
  PatternMining.cpp
)

add_compile_options(-Wno-sign-compare -fopenmp)
add_link_options(-Wno-sign-compare -fopenmp)
add_executable(patternmining ${sources})
add_dependencies(apps patternmining)
target_link_libraries(patternmining PRIVATE Galois::shmem lonestar)
install(TARGETS patternmining DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT apps EXCLUDE_FROM_ALL)

