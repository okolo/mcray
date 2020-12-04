/*
 * main.cpp
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

