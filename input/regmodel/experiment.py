from netpop.compress.experiment import *
import json

if __name__ == '__main__':
	with open('ecoli_exp.json') as f:
		exp = json.load(f)
	with open('ecoli_genome.json') as f:
		genome = json.load(f)

	sim = RegModelExperiment( exp, genome, 'RegModel' )
	
	sim.start()
	sim.advance(1.0)
	sim.finalize()
	res = sim.results()
	
	print res['performance']
	
	print res['result']
	
	print sim.params()