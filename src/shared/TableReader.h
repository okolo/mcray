/*
 * TableReader.h
 *
 *  Created on: Sep 29, 2011
 *      Author: ok
 */

#ifndef TABLEREADER_H_
#define TABLEREADER_H_

#include <istream>
#include <string>
#include <vector>
#include "Utils.h"

namespace Utils {

template<typename X = double >
class TableReaderX : public SmartReferencedObj {
public:
	TableReaderX(std::string aFile, int nColumns=2);
	TableReaderX(std::istream& aInput, int nColumns=2);
	virtual ~TableReaderX();
	inline std::vector<X>& getColumn(int i){return *(m_columns[i]);};
	inline const std::vector<X>& getColumn(int i) const {return *(m_columns[i]);};
protected:
	void Init(std::istream& aInput, int nColumns);
	std::vector<std::vector<X>* > m_columns;
};

typedef TableReaderX<double> TableReader;

} /* namespace Utils */
#endif /* TABLEREADER_H_ */
