/*
 * TableReader.cpp
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


#include "TableReader.h"
#include <iterator>
#include <fstream>
#include <string.h>

namespace Utils {

	template<typename X> TableReaderX<X>::TableReaderX(std::string aFile, int nColumns) {
		std::ifstream file(aFile.c_str());
		if (file.eof())
			Exception::Throw("Failed to open file " + aFile);
		try{
			init(file, nColumns);
		}
		catch(Exception* ex){
			std::string reason = ex->Message();
			delete ex;
			Exception::Throw("Failed to parse " + aFile + " : " + reason);
		}
	}

	template<typename X> TableReaderX<X>::TableReaderX(std::istream& aInput, int nColumns)
	{
		init(aInput, nColumns);
	}

	template<typename X> void TableReaderX<X>::init(std::istream& aInput, int nColumns)
	{
		std::string line;
		bool firstRow=true;
		while (std::getline(aInput, line))
		{
			std::istringstream iss(line);
			X val;
			size_t col=0;
			for(; iss >> val; col++){
				if(firstRow)
					m_columns.push_back(new std::vector<X>());
				else if(col>=m_columns.size())
					Exception::Throw("variable number of records in a row");
				m_columns[col]->push_back(val);
			}
			if(col!=m_columns.size())
				Exception::Throw("variable number of records in a row");
			if(firstRow && nColumns!=Auto && col!=nColumns){
				Exception::Throw("unexpected number of records: " + ToString(col) +  " in a row (" + ToString(nColumns) + " expected)");
			}
			firstRow=false;
		}

//		for(int i=0; i<nColumns; i++)
//			m_columns.push_back(new std::vector<X>());
//
//		std::istream_iterator<X> eos;
//		std::istream_iterator<X> iit (aInput);
//
//		int recNo=0;
//		for(; iit!=eos; recNo++, iit++)
//			m_columns[recNo%nColumns]->push_back(*iit);
//
//		if(recNo%nColumns)
//			Exception::Throw("unexpected number of rows");
	}

	template<typename X> TableReaderX<X>::~TableReaderX()
	{
		for(int i=m_columns.size()-1; i>=0; i--)
			delete m_columns[i];
	}

	template class TableReaderX<double>;
	template class TableReaderX<long double>;

	CTableWithHeaderReader::CTableWithHeaderReader(const char* theFileName, char theCommentSimbol):
			m_commentSymbol(theCommentSimbol),
			m_curHeaderLine(0),
			m_curDataLine(0)
	{
		m_file = fopen(theFileName,"rt");
	}

	CTableWithHeaderReader::~CTableWithHeaderReader()
	{

	}

	bool CTableWithHeaderReader::readLine()
	{
		for(int i=0;i<TAB_READER_BUF_LENGTH;i++)
		{
			int ch = fgetc(m_file);
			if(ch == EOF)
			{
				if(m_buffer[0] == m_commentSymbol)
					i = 0;
				m_buffer[i] = '\0';
				return i;
			}
			if(ch == '\n')
			{
				if((i==0)||(m_buffer[0] == m_commentSymbol))
				{
					i = -1;
					continue;//skipping empty strings and comments
				}
				else
				{
					m_buffer[i] = '\0';
					return true;
				}
			}
			else
				m_buffer[i] = ch;
		}
		Exception::Throw("CDataReader : buffer overflow");
		return false;
	}

	bool CTableWithHeaderReader::read()
	{
		if(!m_file)
			return false;
		//	bool result = false;
		while(readLine())
		{
			m_curHeaderLine++;
			ECompleteStatus status = readHeaderLine(m_buffer);
			switch(status)
			{
				case notCompletedE: continue;
				case failedE: return false;
				case completedE:;
			}
			break;
		}
		while(readLine())
		{
			m_curDataLine++;
			if(!readDataLine(m_buffer))
				return false;
		}
		if(!testData())
			return false;
		processData();
		fclose(m_file);
		m_file = NULL;
		return true;
	}

/// reads number from standard file header line "data length <number>"
/// here number must indicate number of data lines below
	bool CTableWithHeaderReader::readDataLength(int& length, const char* theString)
	{
		const char format[] = " data length %d";
		return (sscanf(theString,format,&length)==1);
	}


/// parameter vector must not be created previously (using CVector::Create() method)
/// default delimiter string " /t/r/n"
	bool CTableWithHeaderReader::readVector(std::vector<double>& v, const char* theString, const char* delimiterStr)
	{
		bool result = false;
		char* stringCopy = 0;
		try{
			const char defaultDelimStr[] = " \t\n\r";
			if (delimiterStr==NULL)
			{
				delimiterStr = defaultDelimStr;
			}

			int length = strlen(theString);
			stringCopy = new char[length + 1];
			strcpy(stringCopy, theString);

			int i=0;
			char* curStr = strtok(stringCopy,delimiterStr);

			// calculating number of entries
			while(curStr)
			{
				double val=0;
				if (sscanf(curStr,"%lg",&val)!=1)
					Exception::Throw("number format error");
				v.push_back(val);
				curStr = strtok(NULL, delimiterStr);
				i++;
			}
			if(i==0)
			{
				Exception::Throw("no tokens");
			}

			result = true;
		}
		catch(const char* aError)
		{
			//error saved and can be retrieved using LastError() method
		}
		delete[] stringCopy;
		return result;
	}

	double CTableWithHeaderReader::yValue(double xValue, double *xArray, double *yArray, int xArraySize, double leftValue, double rightValue)
//find the nearest to xValue x-array element (increasing x-array supposed)
//like CVector::findX function, and than using 1-st order approximation to find y value
	{
		int left = 0;
		int right = xArraySize-1;

		if(xValue>xArray[right])
			return rightValue;
		if(xValue<xArray[left])
			return leftValue;

		for(int i=(left+right)/2;right-left>1;i=(left+right)/2)
		{
			if(xValue>xArray[i])
				left=i;
			else
				right=i;
		}//finding nearest point

		ASSERT((right - left) == 1);

		double result = yArray[left]+(xValue-xArray[left])*(yArray[right]-yArray[left])/(xArray[right]-xArray[left]);
		return result;
	}

	void CTableWithHeaderReader::inverseArray(double *pArray, unsigned int size)
	{
		unsigned int maxIndex = size/2;
		unsigned int i1,i2;
		for(i1=0,i2=size-1;i1<maxIndex;i1++,i2--)
		{
			double mem = pArray[i1];
			pArray[i1] = pArray[i2];
			pArray[i2] = mem;
		}
	}

} /* namespace Utils */
