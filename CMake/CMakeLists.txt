# The name of our project is "GPPG_PROJ".  CMakeLists files in this project can
# refer to the root source directory of the project as ${GPPG_PROJ_SOURCE_DIR} and
# to the root binary directory of the project as ${GPPG_PROJ_BINARY_DIR}.
project (GPPG_PROJ)

# Recurse into the "GPPG" and "CPGSimulator" subdirectories.  This does not actually
# cause another cmake executable to run.  The same process will walk through
# the project's entire directory structure.
add_subdirectory (GPPG)
add_subdirectory (CPGSimulator)
