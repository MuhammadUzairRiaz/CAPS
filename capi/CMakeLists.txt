add_library(caps_c SHARED src/caps_c.cpp)
target_include_directories(caps_c PUBLIC include)
target_link_libraries(caps_c PRIVATE caps_core)
# One unversioned library (libcaps.dylib, libcaps.so, caps.dll): the Studio loads it by name from its own folder.
set_target_properties(caps_c PROPERTIES OUTPUT_NAME caps WINDOWS_EXPORT_ALL_SYMBOLS ON)
if(APPLE)
  set_target_properties(caps_c PROPERTIES INSTALL_NAME_DIR "@rpath" BUILD_WITH_INSTALL_NAME_DIR ON)
endif()
