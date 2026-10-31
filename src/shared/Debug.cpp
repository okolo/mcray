/* 
 * File:   Debug.cpp
 * Author: ok
 * 
 * Created on September 20, 2011, 10:36 AM
 */

#include "Debug.h"
#include "Utils.h"
#include <iostream>
#include <gsl/gsl_errno.h>
#include <fstream>

namespace Utils {

static Debug debug;

bool Debug::fInitialized = false;

Debug::Debug()
{
	if(fInitialized) return;
	gsl_set_error_handler (GslFrrorHandler);
	fInitialized = true;
}

std::ofstream Debug::fOut;

void Debug::SetOutputFile(std::string aFile)
{
#pragma omp critical
	{
		if(fOut.is_open())
			fOut.close();
		fOut.open(aFile.c_str(),std::ios::out);
	}
}

void Debug::PrintLine(std::string aMessage)
{
#pragma omp critical
	{
		if(fOut.is_open())
			fOut << aMessage << std::endl;
		//else
		//	std::cerr << aMessage << std::endl;
	}
}

void Debug::GslFrrorHandler (const char * reason,
        const char * file,
        int line,
        int gsl_errno)
{
	std::ostringstream message;
	message << "GSL error " << gsl_errno << " occurred in " << file << "("<< line << ")\nreason:" << reason << "\n";
	std::string msg = message.str();
	DebugBreakpoint();
	Exception::Throw(msg);
}

void Debug::DebugBreakpoint()
{
/// add breakpoint here
    std::cerr << "debug breakpoint reached\n";
}

}//end of namespace Utils
