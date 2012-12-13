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
#ifndef DOUBLE_DENSE_H
#define DOUBLE_DENSE_H

#ifdef __cplusplus

#include <iostream>
#include <stdexcept>
#include <string.h>
#include <MatrixBase.h>


extern "C" void dgetrf_( const int * M,  const int* N, double* A, const int *lda, int* ipiv, int* result );
extern "C" void dgetrs_( const char* TRANS,  const int* N, const int* nrhs, double* A, const int *lda, int* ipiv, double* B, const int* ldb,  int* result, int tlen);

namespace BMatrix{

/*---------------specialization for double Dense matrices (Dense<double>)-------------*/
template<> Dense<double>::Dense(int m, int n){
	this->rows=m;
	this->cols=n;
	
	this->data = new double[m*n];
	bzero(this->data,m*n*sizeof(double));
	
	this->LU_factors = NULL; //we do not need to allocate it unless we are going to invert the matrix
	have_LU_factors = false;
}

template<> 
void Dense<double>::create(int m, int n){ 
	this->rows=m;
	this->cols=n;
	this->data = new double[m*n];
	this->LU_factors = NULL; //we do not need to allocate it unless we are going to invert the matrix
	have_LU_factors = false;
}
	
template<> Dense<double>& Dense<double>::operator = (const Dense<double> &A){
		
		//if it is a just created matrix, then create it with the same size of A
		if(!this->data){
		    create(A.rows,A.cols);
		}
		if(this->rows!=A.rows || this->cols!=A.cols){
		    throw std::runtime_error("Can not equate two matrices with different sizes");
		}
		
		memcpy(this->data , A.data, this->rows*this->cols*sizeof(double));
	
		LU_factors = NULL; //we are going to allocate it if we are going to invert the matrix
		have_LU_factors = false;
		
		return *this;
}  
	
//this function finds (this)^{-1}*RHS using LU factorization and F/B subst.
///Note: RHS is overwritten with the results 
template<> Dense<double> Dense<double>::solve (double* _RHS , const int Nrhs){
	
	int result;
	int* i_piv = new int[this->rows];
	
	double* RHS = new double[this->cols];
	memcpy(RHS , _RHS , this->cols*sizeof(double));
	
	if(!have_LU_factors){
	    LU_factors = new double[this->rows*this->cols];	
	    memcpy(LU_factors,this->data,this->rows*this->cols*sizeof(double));
	    
	    result = 0;
	    dgetrf_ ( &this->rows, &this->cols, LU_factors, &this->rows, i_piv, &result );	
	    if ( result ){
		    throw std::runtime_error("Cannot find LU factors of a dense matrix");
	    }
	 }

	 result = 0;
	 char tran = 'N';
	 dgetrs_(&tran, &this->rows, &Nrhs, LU_factors, &this->rows, i_piv, RHS, &this->rows,&result , 1);
	 
	 if ( result ){
		  throw std::runtime_error("Cannot Forward/Backward subst for a dense matrix");
	 }
	 
	 delete[] i_piv;
	 
	 return Dense<double>(this->rows, Nrhs , RHS);
}

///Note: RHS is overwritten with the results 
template<> Dense<double> Dense<double>::solve (const Dense<double>& _RHS){

	int result;
	int* i_piv = new int[this->rows];

	Dense<double> RHS = _RHS;
	
	if(!have_LU_factors){
	    LU_factors = new double[this->rows*this->cols];
	    memcpy(LU_factors,this->data,this->rows*this->cols*sizeof(double));
	    result = 0;
	    dgetrf_ ( &this->rows, &this->cols, LU_factors, &this->rows, i_piv, &result );	
	    if ( result ){
		    throw std::runtime_error("Cannot find LU factors of a dense matrix");
	    }
	 }

	 result = 0;
	 char tran = 'N';
	 int Nrhs = RHS.cols;
	 dgetrs_(&tran, &this->rows, &Nrhs, LU_factors, &this->rows, i_piv, RHS.data, &this->rows,&result , 1);
	 
	 if ( result ){
		  throw std::runtime_error("Cannot Forward/Backward subst for a dense matrix");
	 }
	 
	 delete[] i_piv;
	 
	 return RHS;
}

template<> DBase<double>*  Dense<double>::solve(const DBase<double> &B){
	const Dense<double>* BB = dynamic_cast<const Dense<double>*>(&B);
    
	if(BB){
	    return new Dense<double>(this->solve( (*BB) ));
	}else{
	    throw std::runtime_error("I can only handle solve of Dense+Dense");
	}
}

//// result = inv(A)*B
template<> Dense<double> Dense<double>::operator / (const Dense<double> &B){
    if(this->cols!=B.rows){
      throw std::runtime_error("Can not divide two matrices with different rows and columns");
    }
    
    //BMatrix::Dense<double> result(B.rows,B.cols);
    
    //result = B;

    return solve(B);
    
}



}

#endif
#endif