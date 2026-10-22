add_library(caps_core STATIC
  src/live.cpp
  src/update.cpp
  src/mlpot.cpp
  src/ffcomplete.cpp
  src/cluster_job.cpp
  src/elements.cpp
  src/labels.cpp
  src/lattice.cpp
  src/system.cpp
  src/tags.cpp
  src/piece.cpp
  src/bondrules.cpp
  src/probe.cpp
  src/io_lammps.cpp
  src/io_gro.cpp
  src/io_pdb.cpp
  src/io_mol2.cpp
  src/io_car.cpp
  src/io_binary.cpp
  src/io_sdf.cpp
  src/io_detect.cpp
  src/bonds.cpp
  src/analysis.cpp
  src/properties.cpp
    src/mechanics.cpp
    src/mech_props.cpp
    src/cooling_scan.cpp
  src/grow.cpp
  src/semicrystal.cpp
  src/thermal.cpp
  src/field.cpp
  src/nbnxm.cpp
  src/nbnxm_avx2.cpp
  src/power.cpp
  src/relax.cpp
  src/lammps_data.cpp
  src/lammps_read_ff.cpp
  src/gromacs.cpp
  src/chipair.cpp
  src/dynamics.cpp
  src/constraints.cpp
  src/superpose.cpp
  src/equilibrate.cpp
  src/pack.cpp
  src/react.cpp
  src/json.cpp
  src/ffdef.cpp
  src/param_import.cpp
  src/typing.cpp
  src/example_typing.cpp
  src/render.cpp
  src/raytrace.cpp
  src/scene_export.cpp
  src/image.cpp
  src/svg.cpp
  src/smiles.cpp
  src/uff.cpp
  src/crystal.cpp
  src/spacegroup.cpp
  src/solvate.cpp
  src/water.cpp
  src/ffio.cpp
  src/import.cpp
  src/provenance.cpp
  src/voids.cpp
  src/voronoi.cpp
  src/cbmc.cpp
  src/dlpoly.cpp
  src/pairmodel.cpp
  src/adsorption.cpp
  src/slab2d.cpp
  src/adsorb_dft.cpp
  src/vasp_set.cpp
  src/phono3py.cpp
  src/vasp_out.cpp
  src/dft_analysis.cpp
  src/dft_commands.cpp
  src/sorption.cpp
  src/functionalize.cpp
  src/layers.cpp
  src/pipeline_outputs.cpp
  src/ffmerge.cpp
  src/manybody.cpp
  src/dpd.cpp
  src/moltemplate_export.cpp
  src/kremer_grest.cpp
  src/react_view.cpp
  src/bond_react.cpp
  src/reaction_scheme.cpp
  src/io_traj.cpp
  src/yaml.cpp
  src/recipe.cpp
  src/colourvision.cpp
  src/query.cpp
  src/charges.cpp
  src/molinfo.cpp
  src/appearance.cpp
  src/trajectory.cpp
  src/torsion.cpp
  src/normal_modes.cpp
  src/conformers.cpp
  src/ibi.cpp
  src/free_energy.cpp
  src/edit.cpp
  src/edit_script.cpp
  src/interactions.cpp
  src/interface.cpp
  src/nano.cpp
  src/qeq.cpp
  src/checks.cpp
  src/pipeline.cpp
  src/bundle.cpp
  src/pipeline_expr.cpp
  src/pipeline_yaml.cpp
  src/embed.cpp
  src/bench.cpp
  src/polymer.cpp
  src/chimd.cpp
  src/entangle.cpp
  src/polystats.cpp
  src/resolution.cpp
  src/beads.cpp
  src/martini3_molecules.cpp
  src/gromacs_top.cpp
  src/amber.cpp
  src/cg_map.cpp
  src/cg_rules.cpp
  src/cg_commands.cpp
  src/cg_bonded.cpp
  src/cg_nonbonded.cpp
  src/cg_build.cpp
  src/cg_analysis.cpp
  src/cg_backmap.cpp
  src/kspace.cpp
  src/fortran_runtime.cpp
  fortran/caps_kspace.f90
)
target_include_directories(caps_core PUBLIC include)
target_link_libraries(caps_core PUBLIC ZLIB::ZLIB Threads::Threads)
# caps ff complete looks for data/ml/caps-mace-v1.model in the source tree last (a build run from where it was built)
set_source_files_properties(src/ffcomplete.cpp PROPERTIES COMPILE_DEFINITIONS "CAPS_SOURCE_ROOT=\"${PROJECT_SOURCE_DIR}\"")
# the machine-learned potential: only src/mlpot.cpp sees libtorch (its headers need C++20, as the rest of CAPS uses)
if(CAPS_WITH_TORCH)
  find_package(Torch REQUIRED)
  target_link_libraries(caps_core PUBLIC ${TORCH_LIBRARIES})
  set_source_files_properties(src/mlpot.cpp PROPERTIES COMPILE_DEFINITIONS CAPS_WITH_TORCH)
endif()
# the commit the build came from (caps --version, the Studio's check of a cluster's caps): -DCAPS_COMMIT=… for a source
# archive, else git; a new commit re-runs the configure step through the files git updates
if(NOT CAPS_COMMIT)
  execute_process(COMMAND git rev-parse --short=10 HEAD WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                  OUTPUT_VARIABLE CAPS_GIT_COMMIT OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET RESULT_VARIABLE CAPS_GIT_RESULT)
  if(CAPS_GIT_RESULT EQUAL 0 AND CAPS_GIT_COMMIT)
    set(CAPS_COMMIT_VALUE ${CAPS_GIT_COMMIT})
    if(NOT CAPS_COMMIT_DATE)
      # the commit's date (strict ISO 8601): the update check compares it with a release's
      execute_process(COMMAND git show -s --format=%cI HEAD WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                      OUTPUT_VARIABLE CAPS_COMMIT_DATE OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    endif()
    foreach(f HEAD logs/HEAD)
      if(EXISTS ${PROJECT_SOURCE_DIR}/.git/${f})
        set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/.git/${f})
      endif()
    endforeach()
  else()
    set(CAPS_COMMIT_VALUE unknown)
  endif()
else()
  set(CAPS_COMMIT_VALUE ${CAPS_COMMIT})
endif()
if(NOT CAPS_COMMIT_DATE)
  set(CAPS_COMMIT_DATE unknown)   # -DCAPS_COMMIT_DATE=2026-10-22T09:00:00Z for a source archive
endif()
target_compile_definitions(caps_core PRIVATE CAPS_VERSION_STRING="${PROJECT_VERSION}" CAPS_COMMIT_STRING="${CAPS_COMMIT_VALUE}"
                                             CAPS_COMMIT_DATE_STRING="${CAPS_COMMIT_DATE}")
if(NOT MSVC)
  target_compile_options(caps_core PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wpedantic>)
endif()
set_target_properties(caps_core PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fortran_modules)
# The eight-lane cluster kernel: its own file built with AVX2 and FMA on x86-64 (chosen at run time on a processor that
# has them, so the rest stays portable); elsewhere the file is an empty stub
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64)$" AND NOT CMAKE_OSX_ARCHITECTURES MATCHES "arm64")
  if(MSVC)
    set_source_files_properties(src/nbnxm_avx2.cpp PROPERTIES COMPILE_OPTIONS "/arch:AVX2")
  else()
    set_source_files_properties(src/nbnxm_avx2.cpp PROPERTIES COMPILE_OPTIONS "-mavx2;-mfma")
  endif()
endif()
