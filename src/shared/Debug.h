/* 
 * File:   Debug.h
 * Author: ok
 *
 * Created on September 20, 2011, 10:36 AM
 */

#ifndef DEBUG_H
#define	DEBUG_H

#include <string>
#include <sstream>

namespace Utils {
	
template<class T> std::string ToString (const T& aValue)
{
	std::ostringstream logStr;
	logStr << aValue;
	return logStr.str();
}

#define NOT_IMPLEMENTED Exception::Throw(Utils::ToString("not implemented: " __FILE__ "(") + Utils::ToString(__LINE__) +  ")");
	
class Debug{
public:
	Debug();
	static void DebugBreakpoint();
	static void SetOutputFile(std::string aFile);
	static void PrintLine(std::string aMessage);

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
private:
	static void GslFrrorHandler (const char * reason,
        const char * file,
        int line,
        int gsl_errno);
	static bool fInitialized;
	static std::ofstream fOut;
};

}//namespace Utils {





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

