# Make sure the compiler can find include files from our Hello library.
include_directories (${GPPG_PROJ_SOURCE_DIR}/GPPG)

# Make sure the linker can find the Hello library once it is built.
link_directories (${GPPG_PROJ_BINARY_DIR}/GPPG)

# Add executable called "CPGSimulator" that is built from the source files
# "demo.cxx" and "demo_b.cxx".  The extensions are automatically found.
add_executable (CPGSimulator main.cpp)

# Link the executable to the Hello library.
target_link_libraries (CPGSimulator GPPG)
