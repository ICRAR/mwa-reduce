#include <iostream>
#include <fstream>

#include "model.h"

using namespace std;

int main(int argc, char* argv[])
{
	if(argc != 3)
		std::cout << "Syntax: model2kvis <modelfile> <kvis annotation file>\n";
	else {
		Model model(argv[1]);
		ofstream str(argv[2]);
		for(Model::const_iterator srcIter=model.begin(); srcIter!=model.end(); ++srcIter)
		{
			const ModelSource& s = *srcIter;
			const ModelComponent& c = s.Peak();
			str << "TEXT " << c.PosRA()*180.0/M_PI << " " << c.PosDec()*180.0/M_PI << " " << s.Name() << "\n";
		}
	}
}
