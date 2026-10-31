//
// Created by ok on 26/03/2016.
//

#ifndef MCPROPAG_STECKEREBL_H
#define MCPROPAG_STECKEREBL_H

#include "Background.h"

namespace Backgrounds {
    using namespace MCPropag;


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

#endif //MCPROPAG_STECKEREBL_H
