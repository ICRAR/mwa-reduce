#include "imageadder.h"

#include <string>
#include <iostream>

int main(int argc, char *argv[])
{
	if(argc < 5)
	{
		std::cerr << "Syntax: addimg <outimage> <outweights> <inpimage1> <inpbeam1> [<inpimage2> <inpbeam2> ...]\n"
			"All images should be fits files. Use -1 as inpweight for unity weight, or \"-c inpbeam-real inpbeam-imag\" for complex beam.\n";
	}
	const char *outImageName = argv[1];
	const char *outWeightName = argv[2];
	ImageAdder adder;
	for(int argi=3; argi + 1 < argc; argi += 2)
	{
		const char *inpImageName = argv[argi];
		const char *inpWeightName = argv[argi+1];
		if(std::string(inpWeightName) == "-c")
		{
			adder.Add(inpImageName, argv[argi+1], argv[argi+2]);
			argi+=2;
		}
		else {
			adder.Add(inpImageName, inpWeightName);
		}
	}
	std::cout << '\n';
	adder.Finish(outImageName, outWeightName);
}
