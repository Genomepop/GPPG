cmake_minimum_required(VERSION 2.6)

# The name of our project is "GPPG".  CMakeLists files in this project can
# refer to the root source directory of the project as ${GPPG_PROJ_SOURCE_DIR} and
# to the root binary directory of the project as ${GPPG_PROJ_BINARY_DIR}.
project (GPPG)

# Recurse into the "GPPGLib" and "CPGSimulator" subdirectories.  This does not actually
# cause another cmake executable to run.  The same process will walk through
# the project's entire directory structure.
add_subdirectory (GPPGLib)
#add_subdirectory (CPGSimulator)
