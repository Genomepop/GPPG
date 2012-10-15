option(USE_UBIGRAPH "Use ubigraph to visualize evolution" OFF)

set( GPPG_HDR
	GPPG.h
	Base/Genotype.h
	Base/GenotypeFactory.h
	Base/GenotypeHeap.h
	Base/Mutator.h
	Base/Recombinator.h
	Base/Simulator.h
	Model/Pathway/Data.h
	Model/Pathway/Operation.h
	Model/Sequence/Data.h
	Model/Sequence/IO.h
	Model/Sequence/Operation.h
	Operation/BaseCompressionPolicy.h
	Operation/CompressionPolicy.h
	Operation/GreedyLoad.h
	Operation/Operation.h
	Operation/OperationHeap.h
	Operation/Simulator.h
	Simulator/EvoSimulator.h
	Util/binomial.h
	Util/Random.h
	Util/Tools.h
	Util/json/autolink.h
	Util/json/config.h
	Util/json/features.h
	Util/json/forwards.h
	Util/json/json.h
	Util/json/json_batchallocator.h
	Util/json/reader.h
	Util/json/value.h
	Util/json/writer.h
)

set( GPPG_SRC 
	Base/Genotype.cpp
	Base/GenotypeFactory.cpp
	Base/GenotypeHeap.cpp
	Base/Mutator.cpp
	Base/Recombinator.cpp
	Base/Simulator.cpp
	Model/Pathway/Data.cpp
	Model/Pathway/Operation.cpp
	Model/Sequence/Data.cpp
	Model/Sequence/IO.cpp
	Model/Sequence/Operation.cpp
	Operation/BaseCompressionPolicy.cpp
	Operation/CompressionPolicy.cpp
	Operation/GreedyLoad.cpp
	Operation/Operation.cpp
	Operation/OperationHeap.cpp
	Operation/Simulator.cpp
	Simulator/EvoSimulator.cpp
	Util/Random.cpp
	Util/Tools.cpp
	Util/json/json_reader.cpp
	Util/json/json_value.cpp
	Util/json/json_writer.cpp
)

set( EXTRA_LIBS )
# --- Ubigraph support ---
if (USE_UBIGRAPH)
	set (EXTRA_LIBS ${EXTRA_LIBS} xml)
	set( GPPG_HDR ${GPPG_HDR} 
		Util/Ubigraph/ubiclient.h
		Operation/UbiOperationGraph.h
	)
	set( GPPG_SRC ${GPPG_SRC}
		Util/Ubigraph/ubiclient.c
		Operation/UbiOperationGraph.cpp
	)
endif (USE_UBIGRAPH)
# -----------------------

#configure_file( ${GPPG_SOURCE_DIR}/GPPG.h.in ${GPPG_BINARY_DIR}/GPPG.h )

include_directories (${GPPG_SOURCE_DIR}/GPPGLib)

#install (TARGETS GPPG DESTINATION lib)
#install (FILES ${GPPG_HDR} DESTINATION include)
#add_library (GPPG ${GPPG_SRC})
add_executable(CPGSimulator ${GPPG_SOURCE_DIR}/CPGSimulator/main.cpp ${GPPG_SRC})
install (TARGETS CPGSimulator DESTINATION bin)
#target_link_libraries (GPPG ${EXTRA_LIBS})


