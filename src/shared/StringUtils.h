#ifndef __STRINGUTILS_H_
#define __STRINGUTILS_H_

#include <string>
#include <vector>
#include <sstream>

class StringUtils
{
public:
	static int split(const std::string& input, const std::string& delimiter, std::vector<std::string>& results, bool includeEmpties = true);
	template<class T> static bool parse( const std::string& s, T& value)
	{
	   std::istringstream i(s);
	   if (!(i >> value))
	     return false;
	   return true;
	}
};

#endif
