include(FetchContent)
FetchContent_Declare(googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.15.2.tar.gz
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(caps_tests test_core.cpp test_relax.cpp test_dynamics.cpp test_equilibrate.cpp test_pack.cpp test_react.cpp test_ffdef.cpp test_typing.cpp test_properties.cpp test_mechanics.cpp test_molecule.cpp test_polymer.cpp test_kspace.cpp test_uff.cpp test_crystal.cpp test_nano.cpp test_qeq.cpp test_spacegroup.cpp test_solvate.cpp test_appearance.cpp test_trajectory.cpp test_torsion.cpp test_normal_modes.cpp test_conformers.cpp test_edit.cpp test_interactions.cpp test_import.cpp test_image.cpp test_provenance.cpp test_polystats.cpp test_chipair.cpp test_cbmc.cpp test_dlpoly.cpp test_adsorption.cpp test_sorption.cpp test_functionalize.cpp test_dpd.cpp test_layers.cpp test_ffmerge.cpp test_amber.cpp test_slab2d.cpp test_dft_sets.cpp test_dft_analysis.cpp test_cg_rules.cpp test_cg_bonded.cpp test_cg_nonbonded.cpp test_cg_build.cpp test_cg_analysis.cpp test_live.cpp test_cluster_job.cpp test_resume.cpp test_nbnxm.cpp test_ffcomplete.cpp test_water_libraries.cpp test_grow_robust.cpp test_polymer_crystals.cpp test_phono3py.cpp test_thermal.cpp test_update.cpp test_exchange.cpp test_topology.cpp test_exchange_moves.cpp)
target_link_libraries(caps_tests PRIVATE caps_core GTest::gtest_main)
target_compile_definitions(caps_tests PRIVATE CAPS_SAMPLES="${PROJECT_SOURCE_DIR}/samples" CAPS_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
if(CAPS_WITH_TORCH)   # the tests that fit with CAPS's model (FfCompleteMl.*)
  target_compile_definitions(caps_tests PRIVATE CAPS_WITH_TORCH)
endif()
include(GoogleTest)
gtest_discover_tests(caps_tests)

# The Python bindings against the built library (skipped without python3)
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_FOUND)
  add_test(NAME python_bindings COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tests/python/test_bindings.py)
  set_tests_properties(python_bindings PROPERTIES ENVIRONMENT
    "PYTHONPATH=${PROJECT_SOURCE_DIR}/data/python;CAPS_LIB=$<TARGET_FILE:caps_c>;CAPS_SAMPLES=${PROJECT_SOURCE_DIR}/samples")
  # the journal look of the matplotlib figures (caps.plotstyle) and the copy of it in the scripts CAPS writes
  add_test(NAME python_plotstyle COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tests/python/test_plotstyle.py)
  set_tests_properties(python_plotstyle PROPERTIES ENVIRONMENT "PYTHONPATH=${PROJECT_SOURCE_DIR}/data/python;CAPS_LIB=$<TARGET_FILE:caps_c>")
  # the command line and CAPS Studio in step: every caps command and option placed in cli/gui_sync.json, the Studio's
  # "Copy command" options taken by caps, and a React run the same through the C API and the command line
  add_test(NAME cli_gui_sync COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tests/test_cli_gui_sync.py $<TARGET_FILE:caps_cli> ${PROJECT_SOURCE_DIR})
  set_tests_properties(cli_gui_sync PROPERTIES ENVIRONMENT "CAPS_LIB=$<TARGET_FILE:caps_c>")
endif()

# caps job on a fake SLURM cluster (stub sbatch / squeue / sacct / scancel / ws_allocate / ws_release; needs bash)
if(NOT WIN32)
  add_test(NAME cli_job_fake_cluster COMMAND bash ${PROJECT_SOURCE_DIR}/tests/test_cli_job.sh $<TARGET_FILE:caps_cli> ${PROJECT_SOURCE_DIR})
  # caps relax without --ff on a data file caps ff apply wrote: the file's own force field, or a refusal
  add_test(NAME cli_relax_data_own_ff COMMAND bash ${PROJECT_SOURCE_DIR}/tests/test_cli_data_ff.sh $<TARGET_FILE:caps_cli> ${PROJECT_SOURCE_DIR})
  # no tracked file names the conversion program the removed libraries came from (skipped outside a git checkout)
  add_test(NAME repo_ff_sources COMMAND bash ${PROJECT_SOURCE_DIR}/tests/test_repo_ff_sources.sh ${PROJECT_SOURCE_DIR})
  set_tests_properties(repo_ff_sources PROPERTIES SKIP_RETURN_CODE 77)
endif()
