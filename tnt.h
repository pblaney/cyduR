#ifndef tnt_h
#define tnt_h

using namespace std;

//#include <tnt/tnt.h>
#include <tnt/vec.h>
#include <tnt/cmat.h>

using namespace TNT;

#include <vector>

typedef vector< Matrix<int> > MatrixGroup;

template < class T = int > 
class GMatrixGroup : public vector< Matrix<T> >
{
};

namespace TNT {

template <class T>
bool operator==(const TNT::Matrix<T> &A, 
		const TNT::Matrix<T> &B)
{
  bool bResult = true;
  Subscript i,j;
  Subscript M = A.num_rows();
  Subscript N = A.num_cols();

  if(M!=B.num_rows() || N!=B.num_cols())
    bResult = false;

  for (i=0; i<M && bResult; i++)
    for (j=0; j<N && bResult; j++)
       bResult = (A[i][j] == B[i][j]);

  return bResult;
}

template <class T>
bool operator<(const TNT::Matrix<T> &A, 
	       const TNT::Matrix<T> &B)
{
  if(A.num_rows() != B.num_rows())
    return (A.num_rows() < B.num_rows());
  if(A.num_cols() != B.num_cols())
    return (A.num_cols() < B.num_cols());
  const T *a_begin=&(A[0][0]), *a_end=a_begin+A.size();
  const T *b_begin=&(B[0][0]), *b_end=b_begin+B.size();
  return lexicographical_compare(a_begin,a_end,b_begin,b_end);
  //return lexicographical_compare(&A[0],&A[A.size()],&B[0],B[B.size()]);
}

template <class T>
bool operator<(const TNT::Vector<T> &A, 
	       const TNT::Vector<T> &B)
{
  //if(A.size() != B.size()) return (A.size() < B.size());
  return lexicographical_compare(A.begin(),A.end(),B.begin(),B.end());
}

template <class T>
bool operator==(const Vector<T> &A, 
		const Vector<T> &B)
{
  bool bResult = true;

  Subscript i, M = A.size();
  if(M!=B.size())
    bResult = false;

  for (i=0; i<M && bResult; i++)
    bResult = (A[i] == B[i]);

  return bResult;
}
}

template <class T>
bool operator!=(const Vector<T> &A, 
		const Vector<T> &B)
{
  return !(A==B);
}

template <class T>
Matrix<T> duplicateLastColumn(const Matrix<T> &A)
{
  int i,j;
  Matrix<T> B(A.num_rows(),A.num_cols()+1);
  
  for(i=0;i<A.num_rows();i++) {
    for(j=0;j<A.num_cols();j++)
      B[i][j] = A[i][j];
    B[i][A.num_cols()] = A[i][A.num_cols()-1];
  }
  return B;
}

// Haddamard product (component product) of two vectors
template <class T, class X>
Vector<T> hprod(Vector<T> a,Vector<X> b)
{
  int i,n = a.size();
  assert(n==b.size());
  Vector<T> r(n);
  for(i=0;i<n;i++)
    r[i] = a[i]*(T)b[i];
  return r;
}

template <class T>
Matrix<int> discretize(const Matrix<T> &A)
{
  int i,j;
  Matrix<int> B(A.num_rows(),A.num_cols());
  
  for(i=0;i<A.num_rows();i++) {
    for(j=0;j<A.num_cols();j++)
      B[i][j] = (int)A[i][j];
  }
  return B;
}

template <class T>
T sum(Vector<T> a)
{
  int i;
  T r = 0;

  for(i=0;i<a.size();i++) {
    r += a[i];
  }  
  return r;
}

template <class T>
Vector<T> sum(Matrix<T> a)
{
  int i,j;
  Vector<T> r(a.num_cols(),(T)0);

  for(i=0;i<a.num_rows();i++) {
    for(j=0;j<a.num_cols();j++) {
      r[j] += a[i][j];
    }
  }  
  return r;
}

template <class T>
Vector<T> concatenate(Vector<T> a, Vector<T> b)
{
  int i;
  Vector<T> r(a.size()+b.size(),(T)0);
  for(i=0;i<a.size();i++) 
    r[i] = a[i];
  for(i=0;i<b.size();i++) 
    r[a.size()+i] = b[i];
  return r;
}

template <class T>
Matrix<T> cut(const Matrix<T> &x, int i1, int i2, int j1, int j2)
{
  int i,j;
  Matrix<T> r((i2-i1)+1,(j2-j1)+1);
  for(i=i1; i<=i2; i++) {
    for(j=j1; j<=j2; j++)
      r[i-i1][j-j1] = x[i][j];
  }
  return r;
}

template <class T>
Matrix<T> rbind(const Matrix<T> &x, const Vector<T> &y)
{
  int i,j;
  Matrix<T> r;
  if(x.size() == 0)
    r.newsize(1,y.size());
  else {
    assert(x.num_cols() == y.size());
    r.newsize(x.num_rows()+1,x.num_cols());
    for(i=0; i<x.num_rows(); i++) {
      for(j=0; j<x.num_cols(); j++) 
	r[i][j] = x[i][j];
    }
  }
  int lastRow = r.num_rows()-1;
  for(i=0; i<y.size(); i++)
    r[lastRow][i] = y[i];
  return r;
}

#endif


