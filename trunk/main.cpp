 
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>

#include<MatrixBase.h>
#include<Sparse.h>
#include <DoubleSparse.h>
#include <DoubleDense.h>
#include <MWrap.h>


int main(int argc, char *argv[])
{

  printf("Hello, world!\n");

  BMatrix::Dense<double> M(3,3);
  BMatrix::Dense<double> N(3,3);
  BMatrix::Dense<double> Z(3,3);

  M = 2;
  std::cout<<"M=2"<<std::endl;
  std::cout<<M<<std::endl;
  
  M.put(0,0,1);
  M.put(0,1,2);
  M.put(0,2,3);
  
  M.put(1,0,4);
  M.put(1,1,5);
  M.put(1,2,6);
  
  M.put(2,0,7);
  M.put(2,1,8);
  M.put(2,2,10);
  
  N = 2;
  
  std::cout<<"M="<<M<<std::endl;
  std::cout<<"N="<<N<<std::endl;
  
  Z = M+N;
  std::cout<<"M+N="<<Z<<std::endl;
  
   Z = M-N;
  std::cout<<"M-N="<<Z<<std::endl;
  
  Z = M/M;
  std::cout<<"M/M="<<Z<<std::endl;
  


  /*----------Test Sparse------------------------*/

  BMatrix::Sparse<double> A(3,3);
  A.put(0,0,1);
  A.put(1,1,5);
  A.put(2,2,2);
  A.put(2,1,3);
  
  BMatrix::Sparse<double> B(3,3);
  B.put(0,0,1);
  B.put(2,0,5);
  B.put(1,1,2);
  B.put(2,2,3);
  
  std::cout<<"Sparse A :"<< std::endl;
  std::cout<<A;
  


  //L = &C;

/*----------Test BLOCK Matrices------------------------*/	

 BMatrix::MWrap<double> W1(&M) , W2(&N) , W3(&Z);
 W3 = W1+W2;
  
/*----------Test solve functions------------------------*/
BMatrix::Sparse< BMatrix::Dense<double> > C(3000,3000);
BMatrix::Sparse< double > C_double(3000,3000);
BMatrix::Dense< BMatrix::Dense<double> > D(3000,1);
BMatrix::Dense< double > D_double(3000,1);

BMatrix::Dense<double> Temp(1,1);
for(int i=0; i<3000; i++){
  Temp.put(0,0,i+1);
  C.put(i,i,Temp); //Diagonal matrix
  C_double.put(i,i,i+1); //Diagonal matrix
  
  //RHS
  D.put(i,0,Temp);
  D_double.put(i,0,i+1);
}


C.solve(D); //this replaces D
C_double.solve(D_double); //this replaces D

return EXIT_SUCCESS;
}
