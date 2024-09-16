/*
 * Utils.cpp
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


#include "Utils.h"
#include <omp.h>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/dll.hpp>
#include <sstream>
namespace fs = boost::filesystem;

namespace Utils {

    const char PathSeparator =
#ifdef _WIN32
                            '\\';
#else
                            '/';
#endif


    std::string shared_tables_path(){
        std::string tables_folder_name("tables");
        fs::path prog_path = boost::dll::program_location();

        // for some reason boost::filesystem::path.parent_path() and fs::canonical didn't work properly os my Mac OS so we do it manually

        std::stringstream path(prog_path.string());
        std::string segment;
        std::vector<std::string> seglist;

        while(std::getline(path, segment, PathSeparator))
        {
            seglist.push_back(segment);
        }
        
        if(seglist.size() >= 3){
            std::string app_name = seglist.back();
            seglist.pop_back();
            seglist.pop_back();
            seglist.push_back(std::string("shared"));
            seglist.push_back(app_name);
            seglist.push_back(tables_folder_name);
            std::ostringstream os;
            for(std::string& seg : seglist){
                os << seg << PathSeparator;
            }
            std::string table_dir = os.str();
            table_dir.pop_back();
            std::cout << "extracted path:" << table_dir << "\n\n";
            fs::path p = fs::path(table_dir);
            if(fs::exists(fs::status(p)))
                return table_dir;
            else
                std::cout << "type " << fs::status(p).type() << "\n\n";
        }
        return tables_folder_name + PathSeparator;
    }

    int omp_thread_count() {
        int n = 0;
#pragma omp parallel reduction(+:n)
        n += 1;
        return n;
    }

} // end of namespace Utils
