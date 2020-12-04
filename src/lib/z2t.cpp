/*
 * z2t.cpp
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

#include "Cosmology.h"
#include "Units.h"
#include <iostream>
using namespace mcray;

int main(int argc, char** argv) {
    //Test::PrecisionTests::SetPrecision(128);
    int result = 0;
    try
    {
        double z=-1.;
        if(argc<2 || sscanf(argv[1],"%lg",&z)!=1 || z<0) {
            std::cerr << "Calculate light travel time (in Mpc) from remote object\n"
            << "usage: " << argv[0] << "<z>" << std::endl;
            return 1;
        }
        if(!cosmology.IsInitialized())
            cosmology.Init(z+10.);
        double t = (cosmology.getAgeOfUniverse()-cosmology.z2t(z));
        std::cout << t/units.Mpc << std::endl;
    }
    catch(Exception* ex)
    {
        std::cerr << ex->Message() << std::endl;
        result = 1;
    }
    return result;
}