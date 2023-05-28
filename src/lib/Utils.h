/*
 * Utils.h
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


#ifndef UTILS_H
#define	UTILS_H

#include <string>
#include <sstream>
#include <vector>
#include <float.h>
#include <iostream>
#include <iomanip>
#include <omp.h>
#include <sys/types.h>
#include "Debug.h"
#include <algorithm>

namespace Utils {

#define DIR_DELIMITER_STR "/"
#define TABLES_DIR "tables" DIR_DELIMITER_STR

class Exception
{
public:
	Exception(std::string aErrorMessage):fMessage(aErrorMessage){}
	std::string Message() const { return fMessage; }
	inline static void Throw(std::string aErrorMessage) { throw new Exception(aErrorMessage); }
private:
	std::string fMessage;
};

typedef std::vector<double> Vector;

template<class type>
class SafePtr
{
public:
	SafePtr(type* _pType=0):pType(_pType){};
	virtual ~SafePtr(){delete pType;};
	SafePtr& operator=(type* _pType){delete pType; pType=_pType; return (*this);};
	bool isNull(){return (pType==0);};

	operator type*(){return pType;};
	operator const type*() const {return pType;};

	//operator type&(){ASSERT(pType);return *pType;};
	//operator const type&() const {ASSERT(pType);return *pType;};

	type& operator*() {ASSERT(pType);return *pType;};
	const type& operator*() const {ASSERT(pType);return *pType;};

	inline type* operator->() const
	{
		ASSERT(pType);
		return pType;
	}
private:
	type* pType;
};

template<class type>
class AutoDeletePtrArray : public std::vector<type*>
{
public:
	virtual ~AutoDeletePtrArray()
	{
		std::remove_if(std::vector<type*>::begin(),std::vector<type*>::end(),deleteAll);
	};
	static bool deleteAll( type * aElement ) { delete aElement; return true; }

	inline void setIndexShift(int aShift)
	{
		iMin = -aShift;
	};

	inline type& operator()(int aIndex){
		return *std::vector<type*>::at(aIndex - iMin);
	}

	inline const type& operator()(int aIndex) const{
		return *std::vector<type*>::at(aIndex - iMin);
	}

	inline type* operator[](int aIndex){
		return std::vector<type*>::at(aIndex - iMin);
	}

	inline const type* operator[](int aIndex) const{
		return std::vector<type*>::at(aIndex - iMin);
	}

	inline void add(type* aElem){
		std::vector<type*>::push_back(aElem);
	}
private:
	int iMin;
};

	class  ISmartReferencedObj
{
public:
	virtual void addRef()=0;
	virtual void releaseRef()=0;
};

class SmartReferencedObj : public virtual ISmartReferencedObj{
public:
	SmartReferencedObj():iRefCount(0){};
	void addRef(){
		iRefCount++;
	}
	void releaseRef(){
		iRefCount--;
		if (iRefCount<=0) {
			delete this;
		}
	}
	inline int RefCount() const{return iRefCount;};
	virtual ~SmartReferencedObj(){};
private:
	int iRefCount;
};

template <class I>
class TSmartReferencedObj : public I // It is assumed that I is interface which extends ISmartReferencedObj
{
public:
	TSmartReferencedObj():iRefCount(0){};
	void addRef(){
		iRefCount++;
	}
	void releaseRef(){
		iRefCount--;
		if (iRefCount<=0) {
			delete this;
		}
	}
	inline int RefCount() const{return iRefCount;};
	virtual ~TSmartReferencedObj(){};
private:
	int iRefCount;
};

template <class T>//class T should have addRef() & releaseRef() methods
class SmartPtr
{
public:
  SmartPtr(T* pointee = 0) : iPointee(pointee){
		if (iPointee) iPointee->addRef();
  };

  SmartPtr(const SmartPtr<T>& other) : iPointee(other.iPointee){
		if (iPointee) iPointee->addRef();
  };

  inline SmartPtr& operator=(T* pointee){
	  if(iPointee==pointee)
		  return *this;
	  if (iPointee) iPointee->releaseRef();
	  iPointee = pointee;
	  if (pointee) iPointee->addRef();
	  return *this;
  }

  inline SmartPtr& operator=(const SmartPtr<T>& other){
	  if (iPointee) iPointee->releaseRef();
	  iPointee = other.iPointee;
	  if (iPointee) iPointee->addRef();
	  return *this;
  }

  ~SmartPtr(){
		if (iPointee) iPointee->releaseRef();
  }

  inline bool operator==(T* pointee) const{
	  return iPointee==pointee;
  }

  inline T& operator*() const
  {
	return *iPointee;
  }

  inline T* operator->() const
  {
	return iPointee;
  }
  
  inline operator T*() const
  {
	return iPointee;
  }

  inline T& operator[](int) const{
	  //! although operator T* is defined, smart pointers should not be used as C-arrays
	  NOT_IMPLEMENTED;
	  return *iPointee;
  }

private:
  T* iPointee;
};

template<typename T = double >
        class IDataStorageX{
public:
    virtual void save_array(const std::vector<T>& data, const char* name) = 0;
    virtual std::vector<T> load_array(const char* name) = 0;
};

typedef IDataStorageX<double> IDataStorage;

template<typename T = double >
    class BinaryDataStorageX : public IDataStorageX<T>{
    public:
        BinaryDataStorageX(const char* folder):fFolder(folder){};
        void save_array(const std::vector<T>& data, const char* name){
            std::ofstream ofs;
            std::string file = fFolder + "/" + name;
            std::cout << "writing to " << file << std::endl;
            ofs.open(file, std::ios::out | std::ios::binary);
            ofs.write((const char*)data.data(), sizeof(T)*data.size());
            ofs.close();
        }
        std::vector<T> load_array(const char* name){
            std::ifstream ifs;
            std::string file = fFolder + "/" + name;
            std::cout << "reading from " << file << std::endl;
            ifs.open(file, std::ios::in | std::ios::binary);
            //get length of file
            ifs.seekg(0, std::ios::end);
            size_t length = ifs.tellg();
            ifs.seekg(0, std::ios::beg);
            if (length%sizeof(T) != 0)
                Exception::Throw("unexpected size of data");
            std::vector<T> result(length/sizeof(T), 0);
            ifs.read((char*)result.data(), length);
            ifs.close();
            return result;
        }
    private:
        std::string fFolder;
    };

typedef BinaryDataStorageX<double> BinaryDataStorage;

int omp_thread_count();

}//end of namespace Utils

#endif	/* UTILS_H */

