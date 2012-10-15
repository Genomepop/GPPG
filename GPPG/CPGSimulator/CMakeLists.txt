# Make sure the compiler can find include files from our GPPG library.
#include_directories (${GPPG_SOURCE_DIR}/GPPGLib)

# Make sure the linker can find the GPPG library once it is built.
#link_directories (${GPPG_BINARY_DIR}/GPPGLib)

# Add executable called "CPGSimulator" that is built from the source files
#add_executable (CPGSimulator main.cpp)

# Link the executable to the Hello library.
#target_link_libraries (CPGSimulator GPPG)
