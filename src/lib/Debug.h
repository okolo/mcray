/*
 * Debug.h
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

#ifndef DEBUG_H
#define	DEBUG_H

#include <string>
#include <sstream>
#include <fstream>

namespace Utils {
	
template<class T> std::string ToString (const T& aValue)
{
	std::ostringstream logStr;
	logStr << aValue;
	return logStr.str();
}

#define NOT_IMPLEMENTED Exception::Throw(Utils::ToString("not implemented: " __FILE__ "(") + Utils::ToString(__LINE__) +  ")");

    enum LogLevel{
        ErrorLL=0,
        WarningLL,
        MessageLL,
        VerboseLL
    };
	
class Debug{
public:
	Debug(LogLevel aLevel = MessageLL);
	static void DebugBreakpoint();
	void SetOutputFile(std::string aFile);
	void PrintLine(std::string aMessage, LogLevel aLogLevel = MessageLL);
    void EnableTimestamp(bool aEnable = true){
        fTimestamp = aEnable;
    }
    //set to negative to disable thread filtering
    void SetThreadFilter(int aFilter=-1){
        fThreadFilter = aFilter;
    }

    void SetLevel(LogLevel aLevel){
        fLevel = aLevel;
    }

	static inline int IsInfinity(double _val) /* returns nonzero value if argument is Infinity */
	{
		if(_val>1.6e308)
			return 1;
		if(_val<(-1.6e308))
			return -1;
		return 0;
	}
	static int IsValid(double _val) /* returns zero if _val<0 or is Infinity */
	{
		if(_val>=0 && _val < 1.6e308)
		{
			return 1;
		}
		return 0;
	}
    inline bool IsLogEnabled(LogLevel aLevel) { return aLevel<=fLevel && fOut.is_open(); }
private:
	static void GslFrrorHandler (const char * reason,
        const char * file,
        int line,
        int gsl_errno);
	static bool fInitialized;
    bool fTimestamp;
    time_t fStartTime;
    int fThreadFilter;
    LogLevel fLevel;
	std::ofstream fOut;
};

extern Debug debug;

}//namespace Utils {

//macro to ensure log output is calculated only if log is enabled
#define LOG_OUTPUT(str,level) if(Utils::debug.IsLogEnabled(level)) { Utils::debug.PrintLine((str), level);};

#define LOG_VERBOSE(str) LOG_OUTPUT(str,VerboseLL)
#define LOG_MESSAGE(str) LOG_OUTPUT(str,MessageLL)
#define LOG_WARNING(str) LOG_OUTPUT(str,WarningLL)
#define LOG_ERROR(str) LOG_OUTPUT(str,ErrorLL)

///////////////////////////////////////////////////////////////////////////
// Preprocessor definitions for diagnostic support
///////////////////////////////////////////////////////////////////////////
#ifdef _DEBUG

#define ASSERT(f) {if(!(f))\
{\
	Utils::Debug::DebugBreakpoint();\
	Utils::Exception::Throw(Utils::ToString(__FILE__ "(") + Utils::ToString(__LINE__) + ") assertion failed: " #f );\
};}
#define VERIFY(f) ASSERT(f)
#define DEBUG_ONLY(f)      (f)
#define ASSERT_VALID_NO(f) ASSERT(Utils::Debug::IsValid(f))
#define VERIFY_VALID_NO(f) ASSERT(Utils::Debug::IsValid(f))

#else   // _DEBUG

#define ASSERT(f)          ((void)0)
#define VERIFY(f)          ((void)(f))
#define DEBUG_ONLY(f)      ((void)0)
#define ASSERT_VALID_NO(f)          ((void)0)
#define VERIFY_VALID_NO(f)          ((void)(f))

#endif  // _DEBUG


#endif	/* DEBUG_H */

