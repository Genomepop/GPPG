#General Purpose Population Genetic Library (GPPGLib)
GPPGLib is a new standalone C++ components off the shelf (COTS) library for use in population genetic simulators.
Currently, this library contains classes useful for creating individual-based forward time evolutionary simulators --- most notably containing genotype compression algorithms.
Users can use a runtime simulator (CPGSimulator), extend it, or build their own genotype representations.
Subsequent version will break out the simulator from the COTS library; however, in its current form it is easier to bundle and build both together. 

##Versions
* 1.0 - Basic classes for individual-based forward time evolutionary simulators and genotype compression algorithms

##Configure
* download this source
* generate build files with `cmake .` (or replace `.` with an external build directory)
* `USE_UBIGRAPH` supports visualization of the operation graph with the Ubigraph Server

##UBIGRAPH support
Ubigraph is a 3D graph visualization server (http://ubietylab.net/ubigraph/).
GPPGLib supports the visualization of the operation graph over xmlrpc calls to a locally running Ubigraph server (127.0.0.1).
If you include Ubigraph support, you will need the C libraries necessary to make xmlrpc calls, which include: -lxmlrpc_client -lxmlrpc -lxmlrpc_util 
-lxmlrpc_xmlparse -lxmlrpc_xmltok.  
CMake will automatically link to these libraries --- you just need to make sure they are installed.

##Install
In the build directory, do the usual:
* `make`
* `make install` (you will need to use `sudo` if you use the default `/usr/local/bin` installation site)

##Running the Tool
Some example input files are given in the `/input` folder.
To run the simulator, type `CPGSimulator <input file>`.  All options are stored in the json input file.

###Input File
The input file follows JSON (javascript) syntax.  There are several examples which show all the options.

* individuals - integer, size of population
* generations - integer, number of generations
* scaling - number, scaling factor for simulation input
* steps - integer, provides printout of progress per step.  If performance is recorded, then this is the number of steps in the performance recording.
* compression - dictionary, can be `{"name":"Store-Active"}`, `{"name":"Store-Root"}`, or `{"name":"Greedy-Load", "k":50, "t":5}` where k and t are the Greedy-Load parameters
* genotype - dictionary, can be `{"name": "Sequence", "length":100000 }` or `{"name" : "Pathway", "genes" : 300,"tfs" : 300,"regions": [100,300]}`, where genes is the number of genes, tfs is the number of transcription factors, and regions is the range in promoter size.
* operators - list, this depends on the genotype --- look at the examples for the different supported operations
* output - dictionary, if "performance" is a key, then it provides the output file (csv); if "individuals" is provided, then it designates the location to output all the genetic information of each individual