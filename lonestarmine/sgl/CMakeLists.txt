add_executable(sgl sgl.cpp)
target_link_libraries(sgl PRIVATE Galois::shmem lonestar pangolin)
