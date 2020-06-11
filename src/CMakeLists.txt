add_executable("wave2" wave2.cc driver.hh operator.hh bctype.hh norm.hh)
target_link_dune_default_libraries("wave2")

dune_symlink_to_source_files(FILES wave.ini)


add_custom_command(
        TARGET wave2 POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_SOURCE_DIR}/src/domena1.msh
                ${CMAKE_CURRENT_BINARY_DIR})
