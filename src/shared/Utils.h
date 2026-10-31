/* 
 * File:   Utils.h
 * Author: ok
 *
 * Created on September 20, 2011, 10:21 AM
 */

#ifndef UTILS_H
#define	UTILS_H

#include <string>
#include <sstream>
#include <vector>
#include <float.h>
#include <iostream>
#include <iomanip>
//#include <omp.h>
#include "Debug.h"

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

template <class T>
class TSmartReferencedObj : virtual public ISmartReferencedObj, virtual public T{
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

}//end of namespace Utils

#endif	/* UTILS_H */

