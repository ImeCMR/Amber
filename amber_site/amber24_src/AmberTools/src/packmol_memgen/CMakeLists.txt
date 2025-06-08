project(packmol_memgen)

# buildfile for packmol_memgen

# first, build the three native executables it uses
add_subdirectory(packmol_memgen/lib/memembed)
add_subdirectory(packmol_memgen/lib/packmol)
add_subdirectory(packmol_memgen/lib/ppm3)

# --------------------------------------------------------------------
# now, build the Python package

install_python_library(SCRIPT_ARGS --inside-cmake)
