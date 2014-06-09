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
}