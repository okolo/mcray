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
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/dll.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS
#include <sstream>
namespace fs = boost::filesystem;

namespace Utils {

    const char PathSeparator = DIR_DELIMITER_STR[0];

    std::string default_tables_path(){
        // On conda env installations tables are stored in shared folder
        // For standalone installations tables are stored in bin folder
        std::string tables_folder_name("tables");
        fs::path prog_path = boost::dll::program_location();
        fs::path prog_name = prog_path.stem().string();
        fs::path prog_dir = prog_path.parent_path();
        fs::path shared_tables_dir = prog_dir.parent_path() / "shared" / prog_name / tables_folder_name;
        fs::path exe_tables_dir = prog_dir / tables_folder_name;

        if(fs::is_directory(shared_tables_dir))
            return shared_tables_dir.string();
        if(fs::is_directory(exe_tables_dir))
            return exe_tables_dir.string();

        // we don't check for existence of './tables' and
        // don't throw exception here since on some platforms fs::is_directory doesn't work properly
        std::cerr << "WARNING: Neither " << shared_tables_dir.string() << " nor \n" << exe_tables_dir.string() <<
        " was found.\n You can download latest tables from\n" <<
        "https://github.com/okolo/mcray/tree/main/bin/tables\n";

        return tables_folder_name;
    }

    std::string tables_dir = default_tables_path() + DIR_DELIMITER_STR;

    int omp_thread_count() {
        int n = 0;
#pragma omp parallel reduction(+:n)
        n += 1;
        return n;
    }

} // end of namespace Utils
