#include "io.h"
#include "types.h"
#include <sys/stat.h>
#include <sstream>

/// convert string to real number or integer
template <typename T>
T toNumber(std::string str)
{
     T num;
     std::stringstream ss(str); //turn the string into a stream
     ss >> num; //convert
     return num;
}

namespace io
{

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems)
{
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


std::vector<std::string> split(const std::string &s, char delim)
{
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

void makeDirectory(const std::string folderPath)
{
	std::vector<std::string> x = split(folderPath, '/');
	int n = x.size();
	int i = 0;
	std::string folder = x[i];
	mkdir(folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	i++;
	while(i<n)
	{
		folder = folder + '/' + x[i];
		mkdir(folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		i++;
	}
}

}