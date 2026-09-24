add_executable(caps_cli main.cpp)
target_link_libraries(caps_cli PRIVATE caps_core)
set_target_properties(caps_cli PROPERTIES OUTPUT_NAME caps)
target_compile_definitions(caps_cli PRIVATE CAPS_SOURCE_ROOT="${PROJECT_SOURCE_DIR}")
