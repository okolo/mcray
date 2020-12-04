/*
 * Debug.cpp
 *
 * Author:
 *       Oleg Kalashev
 *
 * Copyright (c) 2020 Institute for Nuclear Research, RAS
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

#include "Debug.h"
#include "Utils.h"
#include <iostream>
#include <gsl/gsl_errno.h>
#include <fstream>
#include <omp.h>

#include <execinfo.h>
#include <errno.h>
#include <cxxabi.h>
#include <stdio.h>
#include <signal.h>
#include <stdlib.h>


namespace Utils {

    static inline void printStackTrace( FILE *out = stderr, unsigned int max_frames = 63 )
    {
        fprintf(out, "stack trace:\n");

        // storage array for stack trace address data
        void* addrlist[max_frames+1];

        // retrieve current stack addresses
        unsigned int addrlen = backtrace( addrlist, sizeof( addrlist ) / sizeof( void* ));

        if ( addrlen == 0 )
        {
            fprintf( out, "  \n" );
            return;
        }

        // resolve addresses into strings containing "filename(function+address)",
        // Actually it will be ## program address function + offset
        // this array must be free()-ed
        char** symbollist = backtrace_symbols( addrlist, addrlen );

        size_t funcnamesize = 1024;
        char funcname[1024];

        // iterate over the returned symbol lines. skip the first, it is the
        // address of this function.
        for ( unsigned int i = 4; i < addrlen; i++ )
        {
            char* begin_name   = NULL;
            char* begin_offset = NULL;
            char* end_offset   = NULL;

            // find parentheses and +address offset surrounding the mangled name
#ifdef DARWIN
            // OSX style stack trace
      for ( char *p = symbollist[i]; *p; ++p )
      {
         if (( *p == '_' ) && ( *(p-1) == ' ' ))
            begin_name = p-1;
         else if ( *p == '+' )
            begin_offset = p-1;
      }

      if ( begin_name && begin_offset && ( begin_name < begin_offset ))
      {
         *begin_name++ = '\0';
         *begin_offset++ = '\0';

         // mangled name is now in [begin_name, begin_offset) and caller
         // offset in [begin_offset, end_offset). now apply
         // __cxa_demangle():
         int status;
         char* ret = abi::__cxa_demangle( begin_name, &funcname[0],
                                          &funcnamesize, &status );
         if ( status == 0 )
         {
            funcname = ret; // use possibly realloc()-ed string
            fprintf( out, "  %-30s %-40s %s\n",
                     symbollist[i], funcname, begin_offset );
         } else {
            // demangling failed. Output function name as a C function with
            // no arguments.
            fprintf( out, "  %-30s %-38s() %s\n",
                     symbollist[i], begin_name, begin_offset );
         }

#else // !DARWIN - but is posix
            // not OSX style
            // ./module(function+0x15c) [0x8048a6d]
            for ( char *p = symbollist[i]; *p; ++p )
            {
                if ( *p == '(' )
                    begin_name = p;
                else if ( *p == '+' )
                    begin_offset = p;
                else if ( *p == ')' && ( begin_offset || begin_name ))
                    end_offset = p;
            }

            if ( begin_name && end_offset && ( begin_name < end_offset ))
            {
                *begin_name++   = '\0';
                *end_offset++   = '\0';
                if ( begin_offset )
                    *begin_offset++ = '\0';

                // mangled name is now in [begin_name, begin_offset) and caller
                // offset in [begin_offset, end_offset). now apply
                // __cxa_demangle():

                int status = 0;
                char* ret = abi::__cxa_demangle( begin_name, funcname,
                                                 &funcnamesize, &status );
                char* fname = begin_name;
                if ( status == 0 )
                    fname = ret;

                if ( begin_offset )
                {
                    fprintf( out, "  %-30s ( %-40s  + %-6s) %s\n",
                             symbollist[i], fname, begin_offset, end_offset );
                } else {
                    fprintf( out, "  %-30s ( %-40s    %-6s) %s\n",
                             symbollist[i], fname, "", end_offset );
                }
#endif  // !DARWIN - but is posix
            } else {
                // couldn't parse the line? print the whole line.
                fprintf(out, "  %-40s\n", symbollist[i]);
            }
        }

        free(symbollist);
    }

    void abortHandler( int signum )
    {
        // associate each signal with a signal name string.
        const char* name = NULL;
        switch( signum )
        {
            case SIGABRT: name = "SIGABRT";  break;
            case SIGSEGV: name = "SIGSEGV";  break;
            case SIGBUS:  name = "SIGBUS";   break;
            case SIGILL:  name = "SIGILL";   break;
            case SIGFPE:  name = "SIGFPE";   break;
        }

        // Notify the user which signal was caught. We use printf, because this is the
        // most basic output function. Once you get a crash, it is possible that more
        // complex output systems like streams and the like may be corrupted. So we
        // make the most basic call possible to the lowest level, most
        // standard print function.
        if ( name )
            fprintf( stderr, "Caught signal %d (%s)\n", signum, name );
        else
            fprintf( stderr, "Caught signal %d\n", signum );

        // Dump a stack trace.
        // This is the function we will be implementing next.
        printStackTrace();

        // If you caught one of the above signals, it is likely you just
        // want to quit your program right now.
        exit( signum );
    }

Debug debug;

bool Debug::fInitialized = false;

Debug::Debug(LogLevel aLevel)
{
    fInitialized = false;
    fTimestamp = false;
    fStartTime = time(0);
    fThreadFilter = -1;
    fLevel = aLevel;

	if(fInitialized) return;
	gsl_set_error_handler (GslFrrorHandler);

    signal( SIGABRT, abortHandler );
    signal( SIGSEGV, abortHandler );
    signal( SIGILL,  abortHandler );
    signal( SIGFPE,  abortHandler );

	fInitialized = true;
}

void Debug::SetOutputFile(std::string aFile)
{
#pragma omp critical (Debug)
	{
		if(fOut.is_open())
			fOut.close();
		fOut.open(aFile.c_str(),std::ios::out);
	}
}

void Debug::PrintLine(std::string aMessage, LogLevel aLevel)
{
    int threadId = omp_get_thread_num();
    if(fThreadFilter>=0 && threadId!=fThreadFilter)
        return;

#pragma omp critical (Debug)
	{
		if(IsLogEnabled(aLevel)) {
                if (fTimestamp) {
                    fOut << (time(NULL) - fStartTime) << "\t";
                }
                fOut << threadId << "\t" << aMessage << std::endl;
        }
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
    LOG_ERROR(msg);
#ifdef _DEBUG
    DebugBreakpoint();
#endif
	Exception::Throw(msg);
}

void Debug::DebugBreakpoint()
{
/// add breakpoint here
    std::cerr << "debug breakpoint reached\n";
}

}//end of namespace Utils
