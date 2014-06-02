#include "system.h"

#include <boost/filesystem.hpp>

#include <cstdlib>
#include <fstream>
#include <iostream>

std::string System::FindPythonFilePath(const std::string& filename)
{
	if(boost::filesystem::exists(filename))
		return filename;
	std::cout << "Searching " << filename << "... " << std::flush;
	boost::filesystem::path tempPath = boost::filesystem::unique_path();
	const std::string tempFilename = tempPath.native();  // optional
	std::string command =
		std::string("echo \"import sys\nfor a in sys.path:\n  print a\"|python>") +
		tempFilename;
	system(command.c_str());
	std::ifstream searchPathsFile(tempFilename.c_str());
	if(searchPathsFile.bad())
		throw std::runtime_error(("Error in findPythonFilePath: system call did not create expected temporary file " + tempFilename).c_str());
	while(searchPathsFile.good())
	{
		std::string prefixPath;
		std::getline(searchPathsFile, prefixPath);
		boost::filesystem::path searchPath(prefixPath);
		searchPath /= filename;
		if(boost::filesystem::exists(searchPath))
		{
			const std::string result = searchPath.native();
			std::cout << result << '\n';
			searchPathsFile.close();
			boost::filesystem::remove(tempPath);
			return result;
		}
	}
	searchPathsFile.close();
	boost::filesystem::remove(tempPath);
	throw std::runtime_error(std::string("Could not find Python file ") + filename);
}
