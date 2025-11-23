#ifndef generate_h
#define generate_h

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <set>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

#include "util.h"
#include "tnt.h"
#include "DiscreteHistogram.h"

typedef map<int,int> svType;

void rng_initialise(int);
void rng_free();

void
reportLogLogDegreeDist(const Matrix <discrete> &C);

Matrix<discrete>
createBarabasiNetwork(int);

Matrix<discrete>
createDirectedBarabasiNetwork(int);

Matrix<discrete>
createErdosRenyiNetwork(int,double);

Matrix<discrete>
createUnsignedErdosRenyiNetwork(int,double);

Matrix<discrete>
createUnsignedSiegalBergmanNetwork(int,double);

void
convertUndirectedNetwork(Matrix<discrete> &C);

double 
averageConnectivity(const Matrix<discrete> &C);

template <class T>
int
countNonZeroEntries(const Matrix<T> &C)
{
  int i, j, sum=0, N=C.num_rows();
  FOR(i,N) {
    FOR(j,N)
      if(C[i][j])
	sum++;
  }
  return sum;
}

template <class T>
vector< pair<int,int> >
extractNonZeroEntries(const Matrix<T> &C)
{
  vector< pair<int,int> > vR;
  int i,j, sum=0, N=C.num_rows();
  FOR(i,N) {
    FOR(j,N)
      if(C[i][j])
	vR.push_back( make_pair(i,j) );
  }
  return vR;
}

template <class T>
vector< pair<int,int> >
extractZeroEntries(const Matrix<T> &C)
{
  vector< pair<int,int> > vR;
  int i,j, sum=0, N=C.num_rows();
  FOR(i,N) {
    FOR(j,N)
      if(!C[i][j])
	vR.push_back( make_pair(i,j) );
  }
  return vR;
}

template <class T>
Vector<T>
flattenMatrix(const Matrix<T> &C)
{
  Vector<T> vR(C.size());
  int i,j,k=0, N=C.num_rows();
  FOR(i,N) {
    FOR(j,N)
      vR[k++] = C[i][j];
  }
  return vR;
}
/*
template <class T>
vector<T>
concatenate_vector(const vector<T> &A, const vector<T> &B)
{
  int i,j=0;
  vector<T> vR(A.size()+B.size());
  FOR(i,A.size()) vR[j++] = A[i];
  FOR(i,B.size()) vR[j++] = B[i];
  return vR;
}
*/
template <class T>
Vector<T>
concatenate_vector(const Vector<T> &A, const Vector<T> &B)
{
  int i,j=0;
  Vector<T> vR(A.size()+B.size());
  FOR(i,A.size()) vR[j++] = A[i];
  FOR(i,B.size()) vR[j++] = B[i];
  return vR;
}

template <class T>
Vector<T> 
convertToTNT(const vector<T> &v)
{
  int i;
  Vector<T> vR(v.size());
  copy(v.begin(),v.end(),vR.begin());
  return vR;
}

#endif
