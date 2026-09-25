include(FetchContent)
FetchContent_Declare(googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.15.2.tar.gz
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(caps_tests test_core.cpp test_relax.cpp test_dynamics.cpp test_equilibrate.cpp test_pack.cpp test_react.cpp test_ffdef.cpp test_typing.cpp test_properties.cpp test_mechanics.cpp test_molecule.cpp test_polymer.cpp test_kspace.cpp test_uff.cpp test_crystal.cpp test_nano.cpp test_qeq.cpp test_spacegroup.cpp test_peptide.cpp test_solvate.cpp test_appearance.cpp test_trajectory.cpp test_torsion.cpp test_edit.cpp test_interactions.cpp)
target_link_libraries(caps_tests PRIVATE caps_core GTest::gtest_main)
target_compile_definitions(caps_tests PRIVATE CAPS_SAMPLES="${PROJECT_SOURCE_DIR}/samples" CAPS_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
include(GoogleTest)
gtest_discover_tests(caps_tests)
