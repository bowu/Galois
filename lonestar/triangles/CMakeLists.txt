add_executable(triangles Triangles.cpp)
add_dependencies(apps triangles)
target_link_libraries(triangles PRIVATE Galois::shmem lonestar)
install(TARGETS triangles DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT apps EXCLUDE_FROM_ALL)
add_test_scale(small-node triangles -algo=nodeiterator "${BASEINPUT}/stanford/communities/DBLP/com-dblp.wgt32.sym.gr.triangles")
add_test_scale(small-edge triangles -algo=edgeiterator "${BASEINPUT}/stanford/communities/DBLP/com-dblp.wgt32.sym.gr.triangles")

add_executable(trianglesAOS Triangles_AOS.cpp)
add_dependencies(apps trianglesAOS)
target_link_libraries(trianglesAOS PRIVATE Galois::shmem lonestar)
install(TARGETS trianglesAOS DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT apps EXCLUDE_FROM_ALL)
