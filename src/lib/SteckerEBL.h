/*
 * SteckerEBL.h
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


#ifndef MCRAY_STECKEREBL_H
#define MCRAY_STECKEREBL_H

#include "Background.h"

namespace Backgrounds {
    using namespace mcray;


// IR/O spectrum from Stecker at al. astro-ph/0510449
//(received from Matt Malkan <malkan@astro.UCLA.EDU> on the 6th of Jan 2006)

    class Stecker2005EBL : public CTableWithHeaderReader, public IBackground
    {
    public:
        Stecker2005EBL();
        virtual ~Stecker2005EBL();


        virtual double n(double aE, double aZ) const;
        virtual double MinE() const;
        virtual double MaxE() const;
        virtual double MinZ() const;
        virtual double MaxZ() const;
        virtual std::string Name() const;

        static int UnitTest();

    private:
// from CDataReader
        ECompleteStatus readHeaderLine(const char* theString);
        bool readDataLine(const char* theString);
        bool testData();
        void processData();

// data
        int m_accuracyE;
        int m_accuracyZ;
        int m_curE;
        std::vector<double> m_zArray;// redshifts
        std::vector<double> m_eArray;// energies
        std::vector<std::vector<double> > m_nArray;// photon concentrations
        std::vector<SafePtr<TableFunction> > m_fE; /// F(E,z_i)
    };
}

#endif //MCRAY_STECKEREBL_H
