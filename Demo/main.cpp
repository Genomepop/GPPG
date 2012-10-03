#include <iostream>
#include <fstream>

using namespace std;

#include "Operation/GreedyLoad.h"
#include <Model/Sequence/Operation.h>
#include <Model/Sequence/IO.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <Base/Simulator.h>
#include <Operation/OperationHeap.h>
#include <Operation/BaseCompressionPolicy.h>
#include <Operation/UbiOperationGraph.h>
#include <Simulator/EvoSimulator.h>

#include <Model/Pathway/Operation.h>
#include <Util/json/json.h>

using namespace GPPG;
using namespace GPPG::Model;
using namespace GPPG::Model::TransReg;
using namespace boost::numeric;

void testSequence() {
	//GreedyLoad gl(3);
	
	ublas::vector<double> distr = ublas::vector<double>(4);
	for (int i=0; i<distr.size(); i++) {
		distr(i) = 1.0/distr.size();
	}
	
	cout << distr << endl;
	SequenceRootFactory factory(10, distr);
	
	SequenceRoot* sr = factory.random();
	SequenceData& sd = *sr->data();
	cout << sd << endl; 
	
	ublas::matrix<double> T(4,4);
    for (unsigned i = 0; i < T.size1(); ++ i)
        for (unsigned j = 0; j < T.size2(); ++ j)
            T (i, j) = 0.25;
	
	SequencePointMutator spm(0.1, T);
	
	cout << spm << endl;
	
	cout << "About to mutate.\n";
	OpSequence* m1 = spm.mutate( *sr );
	cout << "Finished mutate.\n";
	cout << *m1 << endl;
	cout << m1->length() << endl;
	for (int i=0; i<m1->length(); i++) {
		cout << m1->get(i);
	}
	cout << endl;
	cout << "Hello World2!\n";
}

void testSimulator() {
	
	//PopulationSimulator psim( new OperationGraph(new BaseCompressionPolicy(STORE_ACTIVE)) );
	//EvoSimulator psim( new OperationGraph(new BaseCompressionPolicy(STORE_ROOT)) );
	EvoSimulator psim( new OperationGraph(new GreedyLoad(20, 5) ));
	
	ublas::vector<double> distr = ublas::vector<double>(4);
	for (int i=0; i<distr.size(); i++) {
		distr(i) = 1.0/distr.size();
	}
	double scaling = 1e2;
	long N = 1e4;
	long L = 1e6;
	long G = 2e7;
	double u = 1e-9;
	double ud = 1e-9;
	double ui = 1e-9;
	double ur = 1e-9;
	cout << "N = " << N/scaling << endl;
	cout << "G = " << G/scaling << endl;
	cout << "u = " << u*scaling << endl;
	
	SequenceRootFactory factory(L, distr);
	SequenceRoot* sr = factory.random();
	
	ublas::matrix<double> T(4,4);
    for (unsigned i = 0; i < T.size1(); ++ i)
        for (unsigned j = 0; j < T.size2(); ++ j)
            T (i, j) = 0.25;
	
	SequencePointMutator *spm = new SequencePointMutator( u*scaling, T);
	SequenceDeletionMutator *sdm = new SequenceDeletionMutator( ud*scaling , 10, 20);
	SequenceInsertionMutator *sim = new SequenceInsertionMutator( ui*scaling , 10, 20, distr);

	psim.addGenotype( sr, 1.0 );
	psim.addMutator( spm );
	//psim.addMutator(sdm);
	//psim.addMutator(sim);

	
	psim.addRecombinator(new SequenceRecombinator(ur*scaling) );
	int steps = 1000;
	for (int i=0; i<steps; i++) {
		psim.evolve( N/scaling, (G/scaling)/steps );	
		//usleep(50000);
		cout << "Done with " << i << endl;
	}
	//psim.evolve( 500, 10000 );
	
	cout << "Generation: " << psim.clock() << endl;
	set<IGenotype*>::iterator git;
	const set<IGenotype*>& active = psim.activeGenotypes();
	int i=0;
	for (git=active.begin(); git!=active.end(); git++) {
		IGenotype* g = *git;
		cout << "Genotype " << i << "/" << g->order() << ": " << g->frequency() << endl;
		i++;
	}
	
}

void testPathway() {
	EvoSimulator psim( new OperationGraph(new BaseCompressionPolicy(STORE_ROOT)) );
	//EvoSimulator psim( new OperationGraph(new GreedyLoad(20, 10) ));
	
	int numGenes = 100;
	int numTFs = 25;
	int minRegion = 10;
	int maxRegion = 100;
	long N = 1e4;
	long G = 1e6;
	
	double scaling = 1e1;
	double u = 1e-9;
	double gain_rate = u*1e-1;
	double loss_rate = 1e-1;
	int overlap = 5;
	
	GlobalInfo* info = PathwayRootFactory::randomInfo( numGenes, numTFs, minRegion, maxRegion );
	
	PathwayRootFactory factory(*info);
	PathwayRoot* root = factory.random();
	

	psim.addGenotype( root, 1.0 );
	psim.addMutator( new BindingSiteMutator( u*scaling, overlap, std::vector<double>(numTFs, gain_rate*scaling), 
											std::vector<double>(numTFs, loss_rate) ) );
					
	int steps = 1000;
	for (int i=0; i<steps; i++) {
		psim.evolve( N/scaling, (G/scaling)/steps );	
		//usleep(50000);
		cout << "Done with " << i << endl;
	}
	//psim.evolve( 500, 10000 );
	
	cout << "Generation: " << psim.clock() << endl;
	set<IGenotype*>::iterator git;
	const set<IGenotype*>& active = psim.activeGenotypes();
	int i=0;
	for (git=active.begin(); git!=active.end(); git++) {
		IGenotype* g = *git;
		cout << "Genotype " << i << "/" << g->order() << ": " << g->frequency() << endl;
		i++;
	}
}

void runSimulation( EvoSimulator* sim, long N, long G, int steps ) {
	cout << "Running Simulation [N="<<N<<", G="<<G<<"]\n";
	
	for (int i=0; i<steps; i++) {
		sim->evolve( N, G/steps );	
		cout << "Done with " << i << " of " << steps << endl;
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
		
		ublas::vector<double> distr = ublas::vector<double>(abet);
		for (int i=0; i<distr.size(); i++) {
			distr(i) = 1.0/distr.size();
		}
		SequenceRootFactory factory(geno["length"].asInt(), distr);
		g = factory.random();
		
		for (int i=0; i<ops.size(); i++) {
			const Json::Value& gOp = ops[i];
			const string& opName = gOp["name"].asString();
			if( opName == "PointMutation" ) {
				ublas::matrix<double> T(abet, abet);
				for (unsigned ii = 0; ii < T.size1(); ++ ii)
					for (unsigned ij = 0; ij < T.size2(); ++ ij)
						T (ii, ij) = 1.0/abet;
				sim->addMutator( new SequencePointMutator( gOp["rate"].asDouble()*scaling, T) );
			} else if( opName == "Insertion" ) {
				sim->addMutator( new SequenceInsertionMutator( gOp["rate"].asDouble()*scaling, 
															  gOp["size"][(Json::Value::ArrayIndex)0].asInt(), gOp["size"][1].asInt(), distr ));
			} else if( opName == "Deletion" ) {
				sim->addMutator( new SequenceDeletionMutator( gOp["rate"].asDouble()*scaling,
															 gOp["size"][(Json::Value::ArrayIndex)0].asInt(), gOp["size"][1].asInt() ));
			} else if( opName == "Recombination" ) {
				sim->addRecombinator( new SequenceRecombinator( gOp["rate"].asDouble()*scaling ));
			}
		}
	} 
	
	else if (genoName == "Pathway") {
		int numGenes = geno["genes"].asInt();
		int numTFs = geno["tfs"].asInt();
		const Json::Value& regions = geno["regions"];
		int minRegion = regions[(Json::Value::ArrayIndex)0].asInt();
		int maxRegion = regions[(Json::Value::ArrayIndex)1].asInt();
		GlobalInfo* info = PathwayRootFactory::randomInfo( numGenes, numTFs, minRegion, maxRegion );
		
		PathwayRootFactory factory(*info);
		g = factory.random();
		
		
		for (int i=0; i<ops.size(); i++) {
			const Json::Value& gOp = ops[i];
			const string& opName = gOp["name"].asString();
			if( opName == "BindingSiteMutation" ) {
				double u = gOp["lossRate"][(Json::Value::ArrayIndex)0].asDouble();
				sim->addMutator( new BindingSiteMutator( u*scaling, gOp["overlap"].asInt(), std::vector<double>(numTFs, gOp["gainRate"].asDouble()*scaling), 
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
	
	runSimulation(sim, config["individuals"].asInt(), config["generations"].asInt(), config.get("steps", 100).asInt());
}

int main (int argc, char * const argv[])
{
	//testSequence();
	//testSimulator();
	//testPathway();
	
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

