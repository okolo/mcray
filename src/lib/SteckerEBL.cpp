/*
 * SteckerEBL.cpp
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


#include <algorithm>
#include "SteckerEBL.h"
#include "TableBackgrounds.h"

using namespace mcray;

namespace Backgrounds {

#define IRO_DATA_FILE "iro_stecker2005"
    Stecker2005EBL::Stecker2005EBL():
            CTableWithHeaderReader((tables_dir + IRO_DATA_FILE).c_str()),
    m_accuracyE(0),
    m_accuracyZ(0),
    m_curE(0)
{
    read();
}

Stecker2005EBL::~Stecker2005EBL()
{

}

double Stecker2005EBL::n(double aE, double aZ) const
{
    aE /= units.eV;//convert internal units to eV
    int iZ=0;
    //	int iE=-1;
    if(aZ < 1e-4)//z=0
    {
        return m_fE[0]->f(aE);// logarithmic approximation on E
    }
    //iZ=m_zArray.findLeftX(z);
    //std::vector<double>::iterator it=std::find_if(m_zArray.begin(),m_zArray.end(),[](double x){return (x > z);});
    iZ=TableFunction::FindLeftX(m_zArray, aZ);

    if(iZ<0)
        return 0.;

    double n1 = m_fE[iZ]->f(aE);// logarithmic approximation on E
    double n2 = m_fE[iZ+1]->f(aE);
    double z1 = m_zArray[iZ];
    double z2 = m_zArray[iZ+1];

    double result = n1 +  (n2-n1)/(z2-z1)*(aZ - z1);  // simple linear approximation on z

    return result;
}

std::string Stecker2005EBL::Name() const
{
    return "Stecker05";
}

double Stecker2005EBL::MinE() const
{
    return m_eArray[0]*units.eV;
}

double Stecker2005EBL::MaxE() const
{
    return m_eArray[m_accuracyE-1]*units.eV;
}

double Stecker2005EBL::MaxZ() const
{
    return m_zArray[m_accuracyZ-1];
}

double Stecker2005EBL::MinZ() const
{
    return 0;
}

CTableWithHeaderReader::ECompleteStatus Stecker2005EBL::readHeaderLine(const char* theString)
{
    if (m_accuracyE<=0)
    {
        return readDataLength(m_accuracyE,theString)?notCompletedE:failedE;
    }

    if (!readVector(m_zArray, theString))
    {
        return failedE;
    }

    m_eArray.insert(m_eArray.begin(),m_accuracyE,0.0);
    m_accuracyZ = m_zArray.size();
    std::vector<double> v;
    m_nArray.insert(m_nArray.begin(),m_accuracyZ,v);
    for(size_t i=0; i<m_accuracyZ; i++)
        m_nArray[i].insert(m_nArray[i].begin(), m_accuracyE, 0.0);

    SafePtr<TableFunction> ptr;
    m_fE.insert(m_fE.begin(),m_accuracyZ, ptr);

    return completedE;
}

bool Stecker2005EBL::readDataLine(const char* theString)
{
    if (m_curE>=m_accuracyE)
    {
        return false;
    }
    std::vector<double> v;
    if(!readVector(v, theString))
        return false;
    ASSERT(v.size()==m_accuracyZ+1);
    double E = pow(10.,v[0]-9.)/241803.; // first column contains log10(E/Hz), converting to eV
    m_eArray[m_curE] = E;
    for(int i=0; i< m_accuracyZ; i++)
        //m_nArray[i][m_curE] = v[i+1]*1e-6*E*pow(1.+m_zArray[i],-3); // multiply by energy and convert m^-3 to cm^-3, multiplied by (1+z)^{-3} to make it in comoving volume
        m_nArray[i][m_curE] = v[i+1]*1e-6/units.cm3; // convert from m^-3 to internal units in physical volume
    m_curE++;
    return true;
}

bool Stecker2005EBL::testData()
{
    return true;
}

void Stecker2005EBL::processData()
{
    for(int i=0; i<m_accuracyZ; i++)
        m_fE[i] = new LinearFunc(m_eArray, m_nArray[i]);
}

int Stecker2005EBL::UnitTest()
{
    Stecker2005EBL backgr;
    //Kneiske1001EBL backgr;
    //CuttedBackground backgr(new Stecker2005EBL(), 0, 0.5*units.eV);
    BackgroundUtils::UnitTest(backgr);
}

}