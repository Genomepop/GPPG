#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

#include "Operation/GreedyLoad.h"
#include <Operation/GreedyLoadMap.h>
#include <Model/Sequence/Operation.h>
#include <Model/Sequence/IO.h>

#include <Base/Simulator.h>
#include <Operation/OperationHeap.h>
#include <Operation/BaseCompressionPolicy.h>
#include <Operation/UbiOperationGraph.h>
#include <Simulator/EvoSimulator.h>

#include <Model/Pathway/Operation.h>
#include <Util/json/json.h>

#define SUPPORTS_RUSAGE

#ifdef SUPPORTS_RUSAGE
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#endif

using namespace GPPG;
using namespace GPPG::Model;
using namespace GPPG::Model::TransReg;


inline double timeToDbl(const timeval& t) { return t.tv_sec + 1.0*t.tv_usec/1e6; }

void recordUsage(ofstream* out, int step, int generation) {
	rusage stats;
	getrusage( RUSAGE_SELF, &stats );
	(*out) << step << "," << generation << "," << (clock()*1.0/CLOCKS_PER_SEC) << ","
			<< timeToDbl(stats.ru_utime) << "," << timeToDbl(stats.ru_stime) << ","
			<< stats.ru_maxrss << "," << stats.ru_ixrss << "," << stats.ru_idrss << "," << stats.ru_isrss << ","
			<< stats.ru_minflt << "," << stats.ru_majflt << ","
			<< stats.ru_nswap << ","
			<< stats.ru_inblock << "," << stats.ru_oublock << ","
			<< stats.ru_msgsnd << "," << stats.ru_msgrcv << ","
			<< stats.ru_nsignals << ","
			<< stats.ru_nvcsw << "," << stats.ru_nivcsw << endl;
	out->flush();
}

// http://linux.die.net/man/2/getrusage

void runSimulation( EvoSimulator* sim, long N, long G, int steps, ofstream* out ) {
	cout << "Running Simulation [N="<<N<<", G="<<G<<"]\n";

#ifdef SUPPORTS_RUSAGE
	if (out) {
		(*out) << "step,gen,wtime,utime,stime,maxrss,ixrss,idrss,isrss,minflt,majflt,nswap,inblock,oublock,msgsnd,msgrcv,nsignals,nvcsw,nivcsw\n";
	}
#endif
	for (int i=0; i<steps; i++) {
		sim->evolve( N, G/steps );	
		cout << "Done with " << i << " of " << steps << endl;
#ifdef SUPPORTS_RUSAGE
		if (out) {
			recordUsage( out, i, sim->clock() );
		}
#endif
	}
	
}

void outputGenotypes( EvoSimulator* sim ) {
	cout << "Generation: " << sim->clock() << endl;
	set<IGenotype*>::iterator git;
	const set<IGenotype*>& active = sim->activeGenotypes();
	int i=0;
	for (git=active.begin(); git!=active.end(); git++) {
		IGenotype* g = *git;
		cout << "Genotype " << i << "/" << g->order() << ": " << g->frequency() << endl;
		i++;
	}
}

EvoSimulator* createSimulator( const Json::Value& config ) {
	double scaling = config.get("scaling",1).asDouble();
	
	// Set Compression
	ICompressionPolicy* policy=0;
	const Json::Value& compression = config["compression"];
	const string& compName = compression.get("name","Store-Root").asString();
	
	if( compName == "Greedy-Load" ) policy = new GreedyLoad(compression.get("k",20).asInt(), compression.get("t",10).asInt());
	else if( compName == "Store-Root" ) policy = new BaseCompressionPolicy(STORE_ROOT); 
	else if( compName == "Store-Active" ) policy = new BaseCompressionPolicy(STORE_ACTIVE);
	else if( compName == "Store-All" ) policy = new BaseCompressionPolicy(STORE_ALL);
	
	if(!policy) {
		cout << "No valid compression policy provided, got: " << compName << endl;
		return 0;
	}
	
	EvoSimulator* sim = new EvoSimulator( new OperationGraph( policy ) );
	
	// Set Factory & Genotype
	const Json::Value& geno = config["genotype"];
	const Json::Value& ops = config["operators"];
	IGenotype* g = 0;
	const string& genoName = geno["name"].asString();

	if( genoName == "Sequence" ) {
		int abet = geno.get("alphabet",4).asInt();
		
		std::vector<double> distr = std::vector<double>(abet, 1.0/abet);

		SequenceRootFactory factory(geno["length"].asInt(), distr);
		g = factory.random();
		
		for (int i=0; i<ops.size(); i++) {
			const Json::Value& gOp = ops[i];
			const string& opName = gOp["name"].asString();
			if( opName == "PointMutation" ) {
				std::vector<double> T(abet*abet);
				for (unsigned ii = 0; ii < abet; ++ ii)
					for (unsigned ij = 0; ij < abet; ++ ij)
						T[ii*abet+ij] = 1.0/abet;
				sim->addMutator( new SequencePointMutator( gOp.get("cost",1).asDouble(), gOp["rate"].asDouble()*scaling, T) );
			} else if( opName == "Insertion" ) {
				sim->addMutator( new SequenceInsertionMutator(gOp.get("cost",10).asDouble(), gOp["rate"].asDouble()*scaling, 
															  gOp["size"][(Json::Value::ArrayIndex)0].asInt(), gOp["size"][1].asInt(), distr ));
			} else if( opName == "Deletion" ) {
				sim->addMutator( new SequenceDeletionMutator(gOp.get("cost",10).asDouble(), gOp["rate"].asDouble()*scaling,
															 gOp["size"][(Json::Value::ArrayIndex)0].asInt(), gOp["size"][1].asInt() ));
			} else if( opName == "Recombination" ) {
				sim->addRecombinator( new SequenceRecombinator(gOp.get("cost",100).asDouble(), gOp["rate"].asDouble()*scaling ));
			}
		}
	} 
	
	else if (genoName == "Pathway") {
		int numGenes = geno["genes"].asInt();
		int numTFs = geno["tfs"].asInt();
		const Json::Value& regions = geno["regions"];
		int minRegion = regions[(Json::Value::ArrayIndex)0].asInt();
		int maxRegion = regions[1].asInt();
		GlobalInfo* info = PathwayRootFactory::randomInfo( numGenes, numTFs, minRegion, maxRegion );
		
		PathwayRootFactory factory(*info);
		g = factory.random();
		
		
		for (int i=0; i<ops.size(); i++) {
			const Json::Value& gOp = ops[i];
			const string& opName = gOp["name"].asString();
			if( opName == "BindingSiteMutation" ) {
				double u = gOp["lossRate"][(Json::Value::ArrayIndex)0].asDouble();
				sim->addMutator( new BindingSiteMutator( gOp.get("cost",1).asDouble(), u*scaling, gOp["overlap"].asInt(), std::vector<double>(numTFs, gOp["gainRate"].asDouble()*scaling), 
														std::vector<double>(numTFs, gOp["lossRate"][1].asDouble()) ) );
			}
		}
	}
	
	if( !g ) {
		cout << "No valid genotype provided, got: " << genoName << endl;
		return 0;
	}
	sim->addGenotype( g, 1.0 );

	
	// Return the simulator!
	return sim;
}

void createAndRunSimulation( const Json::Value& config ) {
	EvoSimulator* sim = createSimulator( config );

	if (!sim) {
		cout << "Failed to create simulator\n";
		return;
	}
	
	const Json::Value& output = config["output"];
	ofstream* perfFile = 0;
	if( output.isMember("performance") ) {
		perfFile = new ofstream();
		perfFile->open( output["performance"].asCString() );
	}
	
	runSimulation(sim, config["individuals"].asInt(), config["generations"].asInt(), config.get("steps", 100).asInt(), perfFile);
	
	if(perfFile) {
		perfFile->close();
		delete perfFile;
	}
	
	if( output.isMember("population") ) {
		outputGenotypes( sim );
	}
	
	delete sim;
}

int main (int argc, char * const argv[])
{
	// First argument is config file.
	if( argc < 2 ) {
		cout << "Please provide a config file\n";
		return -1;
	}
	
	cout << "Reading configuration file " << argv[1] << endl;
	ifstream t( argv[1] );
	stringstream buffer;
	buffer << t.rdbuf();
	
	Json::Value root;
	Json::Reader reader;
	bool parsingSuccessful = reader.parse(buffer.str(), root);
	if (!parsingSuccessful) {
		cout << "Failed to parse configuration\n" << reader.getFormatedErrorMessages();
		return -1;
	}
	cout << root << endl;
	
	createAndRunSimulation(root);
	
	return 0;
}

