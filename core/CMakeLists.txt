add_library(caps_core STATIC
  src/elements.cpp
  src/system.cpp
  src/io_lammps.cpp
  src/io_gro.cpp
  src/io_pdb.cpp
  src/io_mol2.cpp
  src/io_detect.cpp
  src/bonds.cpp
  src/analysis.cpp
  src/properties.cpp
    src/mechanics.cpp
    src/mech_props.cpp
  src/grow.cpp
  src/field.cpp
  src/relax.cpp
  src/lammps_data.cpp
  src/dynamics.cpp
  src/equilibrate.cpp
  src/pack.cpp
  src/react.cpp
  src/json.cpp
  src/ffdef.cpp
  src/dlfield.cpp
  src/typing.cpp
  src/render.cpp
  src/image.cpp
  src/svg.cpp
  src/smiles.cpp
  src/uff.cpp
  src/crystal.cpp
  src/interface.cpp
  src/nano.cpp
  src/embed.cpp
  src/bench.cpp
  src/polymer.cpp
  src/kspace.cpp
  src/fortran_runtime.cpp
  fortran/caps_kspace.f90
)
target_include_directories(caps_core PUBLIC include)
target_link_libraries(caps_core PUBLIC ZLIB::ZLIB Threads::Threads)
if(NOT MSVC)
  target_compile_options(caps_core PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wpedantic>)
endif()
set_target_properties(caps_core PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fortran_modules)
