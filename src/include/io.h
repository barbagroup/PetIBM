/**
* @file io.h
* @brief Functions for input and output.
*/

#pragma once

#include "types.h"
#include <vector>
#include <string>

/**
* @namespace io
* @brief     Contains all the functions related to I/O
*/
namespace io
{
	std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
	std::vector<std::string> split(const std::string &s, char delim);
	void makeDirectory(const std::string s);

	// output
	//void printSimulationInfo(parameterDB &DB, domain &D);
	//void printTimingInfo(Logger &logger);
	
	//void writeInfoFile(parameterDB &DB, domain &D);
	//void writeGrid(std::string &caseFolder, domain &D);
	//template <typename Vector>
	//void writeData(std::string &caseFolder, int n, Vector &q, Vector &lambda, domain &D);//, bodies &B);
	//void printDeviceMemoryUsage(char *label);
}