/*
 * TableReader.cpp
 *
 *  Created on: Sep 29, 2011
 *      Author: ok
 */

#include "TableReader.h"
#include <iterator>
#include <fstream>

namespace Utils {

template<typename X> TableReaderX<X>::TableReaderX(std::string aFile, int nColumns)
{
	std::ifstream file;
	file.open(aFile.c_str(), std::ios::in);
	if(!file.is_open())
		Exception::Throw("Failed to open " + aFile);
	Init(file, nColumns);
}

template<typename X> TableReaderX<X>::TableReaderX(std::istream& aInput, int nColumns)
{
	Init(aInput, nColumns);
}

template<typename X> void TableReaderX<X>::Init(std::istream& aInput, int nColumns)
{
	for(int i=0; i<nColumns; i++)
		m_columns.push_back(new std::vector<X>());

	std::istream_iterator<X> eos;
	std::istream_iterator<X> iit (aInput);

	int recNo=0;
	for(; iit!=eos; recNo++, iit++)
		m_columns[recNo%nColumns]->push_back(*iit);

	if(recNo%nColumns)
		Exception::Throw("unexpected number of rows");
}

template<typename X> TableReaderX<X>::~TableReaderX()
{
	for(int i=m_columns.size()-1; i>=0; i--)
		delete m_columns[i];
}

template class TableReaderX<double>;
template class TableReaderX<long double>;

} /* namespace Utils */
