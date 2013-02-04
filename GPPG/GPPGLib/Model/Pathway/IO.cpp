#include "IO.h"

using namespace std;

ostream& operator<<(ostream& output, const GPPG::Model::TransReg::GlobalInfo& info) {
	output << "Global Info (Genes=" << info.numGenes() << ", Motifs=" << info.numMotifs() << ", TFs=" << info.numTFs() << ")" << endl;
	output << "Binding: " << endl;
	for(int i=0; i<info.numMotifs(); i++) {
		output << info.getMotifName(i) << ": " << info.getMotifSequence( info.getMotifName(i) ) << " [";
		vector<int> binding = info.binding(i);
		for(int j=0; j<binding.size(); j++) {
			output << binding[j] << ", ";
		}
		output << "]" << endl;
	}
	return output;
}