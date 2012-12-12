//
// C++ Interface: MatrixBase
//
// Description: 
//
//
// Author:  <>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MATRIX_BASE_H
#define MATRIX_BASE_H

#ifdef __cplusplus

#include <iostream>
#include <stdexcept>
#include <string.h>

extern "C" void dgetrf_( const int * M,  const int* N, double* A, const int *lda, int* ipiv, int* result );
extern "C" void dgetrs_( const char* TRANS,  const int* N, const int* nrhs, double* A, const int *lda, int* ipiv, double* B, const int* ldb,  int* result, int tlen);

namespace BMatrix{

template <class T>
class DBase{

protected:
	int rows; //number of rows
	int cols; //number of columns		
	T* data;
public:
	int get_number_of_rows(){ return rows;}
	int get_number_of_cols(){return cols;}

	virtual void put(int m, int n, T value)=0; 
	virtual T get(int m, int n) const=0;

	virtual double det()=0; 
	virtual bool is_zero()=0;

	virtual DBase<T>& operator=(T val)=0;
	virtual DBase<T>& operator/=(T val)=0;

	virtual DBase<T>* clone()=0; //clone and return pointer to the new cloned matrix

	virtual DBase<T>& operator +=(const DBase<T> &A)=0;

	virtual DBase<T>& operator -=(const DBase<T> &A)=0;

	virtual DBase<T>& operator *=(const DBase<T> &A)=0;

	//virtual DBase<T>& operator /=(const DBase<T> &A)=0;

	virtual DBase<T>* operator+(const DBase<T> &B) const = 0;
	virtual DBase<T>* operator-(const DBase<T> &B) const = 0;
	virtual DBase<T>* operator*(const DBase<T> &B) const = 0;
	
	virtual DBase<T>* scale (T* val) const = 0;
	virtual DBase<T>* solve(const DBase<T> &B) const = 0;
	
	~DBase(){
	    if(data) delete[] data;
	}
};

template <class T>
class Dense: public DBase<T>{
	
private:  
	T* LU_factors;
	bool have_LU_factors;
public:

	Dense(){ 
		this->rows=0;
		this->cols=0;
		this->data = NULL;
		this->LU_factors = NULL;
		have_LU_factors = false;
		
	}

	Dense(int m, int n){ 
		this->rows=m;
		this->cols=n;
		this->data = new T[m*n];
		this->LU_factors = new T[m*n];
		
		have_LU_factors = false;
	}

	Dense(int m, int n, T* _data){ 
		this->rows=m;
		this->cols=n;
		this->data = _data;
		this->LU_factors = new T[m*n];
		
		have_LU_factors = false;
	}

 	Dense(const Dense<T> &B){ 
		this->rows= B.rows;
		this->cols= B.cols;
		this->data = new T[this->rows*this->cols];
		
		for(int j=0; j<this->cols; j++){
		    for(int i=0; i<this->rows; i++){
			this->data[i+this->rows*j] = B.data[i+this->rows*j];
		    }  
		}
		
		this->LU_factors = new T[this->rows*this->cols];
		have_LU_factors = false;
	}

	Dense<T>* clone(){
		//return new Dense<T>(*this);
		return new Dense(this->rows, this->cols, this->data);
	}

	~Dense(){
	    if(LU_factors) delete[] LU_factors;
	}
	
	void create(int m, int n){ 
		this->rows=m;
		this->cols=n;
		this->data = new T[m*n];
		this->LU_factors = new T[m*n];
		
		have_LU_factors = false;
	}
	
	Dense<T>& operator = (const Dense<T> &A){
		if(this->rows!=A.rows || this->cols!=A.cols){
		    throw std::runtime_error("Can not equate two matrices with different sizes");
		}
		
		for(int j=0; j<this->cols; j++){
		    for(int i=0; i<this->rows; i++){
			this->data[i+this->rows*j] = A.data[i+this->rows*j];
			
		    }  
		}

		this->LU_factors = new T[this->rows*this->cols];
		have_LU_factors = false;
		
		return *this;
	}  
	
	//this operator will be used to fill the whole matrix with one value
	Dense<T>& operator= (T val){
		for(int j=0; j< this->cols; j++){
		    for(int i=0; i< this->rows; i++){
			this->data[i+this->rows*j] = val;
		    }  
		} 
		
		return *this;
	} 

	Dense<T>& operator/= (T val){
		for(int j=0; j< this->cols; j++){
		    for(int i=0; i< this->rows; i++){
			this->data[i+this->rows*j] /= val;
		    }  
		} 
		
		return *this;
	} 
	
	Dense<T>& operator/= (const Dense<T> &A){
		*this = solve(A);
		
		return *this;
	} 
	
	
	
	
	
	Dense<T> operator /(T val){
		Dense<T> result(this->rows, this->cols);
		for(int j=0; j< this->cols; j++){
		    for(int i=0; i< this->rows; i++){
			result.data[i+this->rows*j] = this->data[i+this->rows*j] / val;
		    }  
		} 
  		return result;
	}
	
	DBase<T>* scale (T* val) const{
		Dense<T>* result = new Dense<T> (this->rows, this->cols);
		for(int j=0; j< this->cols; j++){
		    for(int i=0; i< this->rows; i++){
			result->data[i+this->rows*j] = this->data[i+this->rows*j] / (*val);
		    }  
		} 
  		return result;
	}
	
	T get(int m, int n) const{
		return this->data[m+this->rows*n];
	}
	
	T* get_pointer(int m, int n) const{
		return &(this->data[m+this->rows*n]);
	}

	void put(int m, int n, T value){
		this->data[m+this->rows*n] = value;
		
		have_LU_factors = false;
	}

	double det(){
	      
	}

	bool is_zero(){
		//we only check the diagonals
		for(int i=0; i< this->rows; i++){
		    if(this->data[i+this->rows*i]!=0){
			return false;
		    }
		}
		return true;
	}

	Dense<T> solve(T* RHS, const int Nrhs=1) const{
		throw std::runtime_error("Sorry: Can only solve dense<double> matrices now.");
	}
    
	Dense<T> solve(const Dense<T>& _RHS) const{
		throw std::runtime_error("Sorry: Can only solve dense<double> matrices now.");
	}
	
	DBase<T>* solve(const DBase<T> &B) const{
		throw std::runtime_error("Sorry: Can only solve dense<double> matrices now.");
	}
	
	Dense<T> subtract(const Dense<T> &B) const{
		if(this->rows!=B.rows || this->cols!=B.cols){
			throw std::runtime_error("Can not subtract two matrices with wrong different size");
		}

		int m = this->rows;
		int n = this->cols;
		T value;
		Dense<T> results(m , n);
		for(int j=0; j<n; j++){	
			for(int i=0; i<m; i++){
				value = this->data[i+this->rows*j]-B.data[i+B.rows*j];
				results.put(i,j,value);
			}
		}

		return results;
	}

	Dense<T> add(const Dense<T> &B) const{
		if(this->rows!=B.rows || this->cols!=B.cols){
			throw std::runtime_error("Can not add two matrices with different sizes");
		}

		int m = this->rows;
		int n = this->cols;
		T value;
		Dense<T> results(m , n);
		for(int i=0; i<n; i++){
			for(int j=0; j<m; j++){	
				//value = A.data[i+A.rows*j] + B.data[i+B.rows*j];
				value = get(i,j) + B.data[i+B.rows*j];
				results.put(i,j,value);
			}
		}

		return results;
	}

	Dense<T> multibly(const Dense<T> &B) const{
		if (this->cols != B.rows) {
		    throw std::runtime_error("Can not multibly two matrices with different sizes");
		}

		Dense<T> result = Dense<T>(this->rows, B.cols);

		T v;
		for (int i = 0; i < this->rows; i++) {
		    for (int j = 0; j < B.cols; j++) {
			v = 0.0;

			for (int k = 0; k < B.rows; k++) {
			    v += this->get(i, k) * B.get(k, j);
			}

			result.put(i, j, v);
		    }
		}

		return result;	
	}

	//derefernce operator (*A)
	virtual T* operator* () const {return this->data;}
	
	Dense<T>& operator +=(const Dense<T> &A){
		if(this->rows!=A.rows){
    			throw std::runtime_error("Can not add two matrices with different sizes");
  		}
  
  		for(int i=0; i<this->cols; i++){
			for(int j=0; j<this->rows; j++){	
				this->data[i+this->rows*j] += A.data[i+this->rows*j];
			}
		}

  		have_LU_factors = false;
  		return *this;
	}	

	DBase<T>& operator +=(const DBase<T> &A){
  
  		if (const Dense<T>* BB = dynamic_cast<const Dense<T>*>(&A)) {
			if(this->rows!=BB->rows){
    				throw std::runtime_error("Can not add two matrices with different sizes");
  			}
    			for(int i=0; i<this->cols; i++){
				for(int j=0; j<this->rows; j++){	
					this->data[i+this->rows*j] += BB->data[i+this->rows*j];
				}
			}
  		}else{
    			throw std::runtime_error("I can only handle addition of Dense+Dense");
  		}

  		

  		have_LU_factors = false;
  		return *this;
	}		

	Dense<T>& operator -=(const Dense<T> &A){
		if(this->rows!=A.rows){
    			throw std::runtime_error("Can not subtract two matrices with different sizes");
  		}
  
  		for(int i=0; i<this->cols; i++){
			for(int j=0; j<this->rows; j++){	
				this->data[i+this->rows*j] -= A.data[i+this->rows*j];
			}
		}

  		have_LU_factors = false;
  
  		return *this;
	}

	DBase<T>& operator -=(const DBase<T> &A){

  		if (const Dense<T>* BB = dynamic_cast<const Dense<T>*>(&A)) {
			if(this->rows!=BB->rows){
    				throw std::runtime_error("Can not add two matrices with different sizes");
		  	}

    			for(int i=0; i<this->cols; i++){
				for(int j=0; j<this->rows; j++){	
					this->data[i+this->rows*j] -= BB->data[i+this->rows*j];
				}
			}
  		}else{
    			throw std::runtime_error("I can only handle subtraction of Dense+Dense");
  		}

  		

  		have_LU_factors = false;
  		return *this;
	}	
	
	Dense<T>& operator *=(const Dense<T> &B){
		if (this->cols != B.rows) {
		    throw std::runtime_error("Can not multibly two matrices with different sizes");
		}

		for (int i = 0; i < this->rows; i++) {
		    for (int j = 0; j < B.cols; j++) {
			double v = 0.0;

			for (int k = 0; k < B.rows; k++) {
			    v += this->get(i, k) * B.get(k, j);
			}

			this->put(i, j, v);
		    }
		}

		return *this;	
	}
	
	

	DBase<T>& operator *=(const DBase<T> &A){

  		if (const Dense<T>* BB = dynamic_cast<const Dense<T>*>(&A)) {
    			if (this->cols != BB->rows) {
			    throw std::runtime_error("Can not multibly two matrices with different sizes");
			}

			T v;
			for (int i = 0; i < this->rows; i++) {
			    for (int j = 0; j < BB->cols; j++) {
				v = 0.0;

				for (int k = 0; k < BB->rows; k++) {
				    v += this->get(i, k) * BB->get(k, j);
				}

				this->put(i, j, v);
			    }
			}
  		}else{
    			throw std::runtime_error("I can only handle addition of Dense+Dense");
  		}

  		return *this;
	}	

	/*Dense<T>& operator /=(const Dense<T> &A){
		
  		return *this;
	}*/

	Dense<T> operator +(const Dense<T> &B) const{
		return add(B);
	}

	DBase<T>* operator+(DBase<T> const &B) const{
	  
  		if (const Dense<T>* BB = dynamic_cast<const Dense<T>*>(&B)) {
    			return new Dense<T>( *this + *BB  );
  		}else{
    			throw std::runtime_error("I can only handle addition of Dense+Dense");
  		}
	}

	Dense<T> operator -(const Dense<T> &B) const{
		return subtract(B);
	}

	DBase<T>* operator-(DBase<T> const &B) const{
  		const Dense<T>* BB = dynamic_cast<const Dense<T>*>(&B);

  		if (BB) {
    			return new Dense<T>( this->subtract( (*BB) ) );
  		}else{
    			throw std::runtime_error("I can only handle subtraction of Dense+Dense");
  		}
	}

	Dense<T> operator *(const Dense<T> &B) const{
		return multibly(B);
	}


	DBase<T>* operator*(const DBase<T> &B) const{
  		const Dense<T>* BB = dynamic_cast<const Dense<T>*>(&B);

  		if (BB) {
    			return new Dense<T>( this->multibly( (*BB) ) );
  		}else{
    			throw std::runtime_error("I can only handle addition of Dense+Dense");
  		}
	}

	//// this operator gives result = inv(A)*B
	Dense<T> operator /(const Dense<T> &B) const{
		throw std::runtime_error("TODO: division on general dense matrices");
	}
	  
	bool operator != (T val){
		for(int i=0; i< this->rows; i++){
		    for(int j=0; j< this->rows; j++){
			  if(this->data[i+this->rows*j]!=0){
			      return false;
			  }
		    }
		}
		return true;
	}
	template<class U> friend std::ostream& operator<< (std::ostream &out, const Dense<U> &B);
	
};


template<class U> std::ostream& operator<< (std::ostream &out, const Dense<U> &B){
      out<<"["<<std::endl;
      for(int i=0; i<B.rows; i++){
	  for(int j=0; j<B.cols; j++){
	      out<<B.data[i+B.rows*j]<<" ";
	  }
	  out<<std::endl;
      }
      out<<"]";
      return out;
}



/*template <class U> 
Dense<U> operator +(const Dense<U>& A, const Dense<U>& B)
{
	if(A.rows!=B.rows || A.cols!=B.cols){
		throw std::runtime_error("Can not add two matrices with different sizes");
	}

	int m = A.rows;
	int n = A.cols;
	U value;
	Dense<U> results(m , n);
	for(int i=0; i<n; i++){
		for(int j=0; j<m; j++){	
			//value = A.data[i+A.rows*j] + B.data[i+B.rows*j];
			value = A.get(i,j) + B.data[i+B.rows*j];
			results.put(i,j,value);
		}
	}

	return results;
}*/


	
};
#endif
#endif