add_executable(kcl kcl.cpp)
target_link_libraries(kcl PRIVATE Galois::shmem lonestar pangolin)
