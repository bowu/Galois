set(sources
  src/equivalence.cpp
  src/quick_pattern.cpp
  src/base_embedding.cpp
  src/vertex_embedding.cpp
  src/BfsMining/embedding_list.cpp
)

add_library(pangolin STATIC ${sources})

target_link_libraries(pangolin galois_shmem llvm)

target_include_directories(pangolin PUBLIC
  ${CMAKE_SOURCE_DIR}/libllvm/include
  ${CMAKE_SOURCE_DIR}/libgalois/include
  ${CMAKE_SOURCE_DIR}/libpangolin/include
  ${CMAKE_SOURCE_DIR}/lonestar/include/pangolin
  ${CMAKE_CURRENT_SOURCE_DIR}/include
)

set_target_properties (pangolin PROPERTIES
  INTERFACE_POSITION_INDEPENDENT_CODE On
  POSITION_INDEPENDENT_CODE On
)
