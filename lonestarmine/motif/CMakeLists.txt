add_executable(motif motif.cpp)
target_link_libraries(motif PRIVATE Galois::shmem lonestar pangolin)
