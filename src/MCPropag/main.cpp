/* 
 * File:   main.cpp
 * Author: ok
 *
 * Created on August 26, 2011, 3:33 PM
 */

#include <cstdlib>
#include <time.h>
#include "Test.h"

using namespace Utils;

/*
 * user_main function performing user defined tasks
 * the function should be provided by end user
 */
extern int user_main(int argc, char** argv);

int main(int argc, char** argv) {
	int result = 0;
	time_t startTime = time(0);
	try
	{
		if(argc>=2 && (strcmp(argv[1],"--test"))==0) {
			Test::TestList testList;
			result = testList.Main(argc - 2, argv + 2);
		}
		else
			result = user_main(argc, argv);
	}
	catch(Exception* ex)
	{
		std::cerr << ex->Message() << std::endl;
		result = 1;
	}
	if(!result){
		std::cerr << "# calculation took " << time(0)-startTime << " sec" << std::endl;
	}
	return result;
}

