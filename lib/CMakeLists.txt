add_library(common-lib dummy.f90)

target_link_libraries(common-lib)

install(TARGETS common-lib
        EXPORT common-lib_targets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})