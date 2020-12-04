/*
 * TableReader.h
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

#ifndef TABLEREADER_H_
#define TABLEREADER_H_

#include <istream>
#include <vector>
#include "Utils.h"

namespace Utils {

template<typename X = double >
class TableReaderX : public SmartReferencedObj {
public:
	enum{
		Auto=-1
	};
    //use TableReader::Auto for second argument if number of columns is not known in advance
	TableReaderX(std::string aFile, int nColumns /* =Auto */);
	TableReaderX(std::istream& aInput, int nColumns /* =Auto */);
    virtual ~TableReaderX();
	inline std::vector<X>& getColumn(int i){return *(m_columns[i]);};
	inline const std::vector<X>& getColumn(int i) const {return *(m_columns[i]);};
	inline size_t numberOfColumns() const { return m_columns.size(); };
protected:
	void init(std::istream& aInput, int nColumns);
	std::vector<std::vector<X>* > m_columns;
};

typedef TableReaderX<double> TableReader;

#ifndef TAB_READER_BUF_LENGTH
#define TAB_READER_BUF_LENGTH 65536
#endif

	class CTableWithHeaderReader
	{
	public:
		static void inverseArray(double* pArray, unsigned int size);
		static double yValue(double xValue, double* xArray, double* yArray,int xArraySize, double leftValue=0, double rightValue=0);
		enum ECompleteStatus
		{
			notCompletedE = -1,
			failedE = 0,
			completedE = 1
		};
		CTableWithHeaderReader(const char* theFileName, char theCommentSimbol = '#');
		virtual bool read();
		virtual ~CTableWithHeaderReader();

		/// reads number from standard file header line "data length <number>"
		/// here number must indicate number of data lines below
		/// usually called from readHeaderLine
		/// parameter theString is header string passed to readHeaderLine
		static bool readDataLength(int& length, const char* theString);

		/// the data read will be appended to vector
		/// default delimiter string " /t/r/n"
		static bool readVector(std::vector<double>& aVector, const char* aString, const char* aDelimiterStr=NULL);

	protected:
		virtual ECompleteStatus readHeaderLine(const char* theString){return completedE;};

		//readDataLine returns false in case of error
		virtual bool readDataLine(const char* theString) = 0;
		virtual bool testData(){return true;};
		virtual void processData(){};

		char m_commentSymbol;
		int m_curHeaderLine;//starting from 1
		int m_curDataLine;//starting from 1
	private:
		bool readLine();//returns false if EOF was reached

		char m_buffer[TAB_READER_BUF_LENGTH];
		FILE* m_file;
	};

} /* namespace Utils */
#endif /* TABLEREADER_H_ */
