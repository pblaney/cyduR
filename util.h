#ifndef util_h
#define util_h

/*
  Small utility functions
*/

#include <fstream>
#include <list>
#include <string>
#include <map>
#include <iterator>
#include <functional>

#include "typedefs.h"
#include "tnt.h"
#include "StateVector.h"
#include <gsl/gsl_vector.h>

void fatalError(char *errMessage);

void fatalError(string errMessage);

class Exception {};

class IllegalArgumentException: public Exception {};

template <class T> void readFile(T &result, const char *name) {
  ifstream inputFile(name);
  if(!inputFile) 
    fatalError("cannot open file "+string(name));
  inputFile >> result;
  inputFile.close();
}

template <class T> void readFile(T &result, const string sName) {
  readFile(result,sName.c_str());
}

template <class T> void writeFile(T &data, const char *name) {
  ofstream outputFile(name);
  outputFile << data;
  outputFile.close();
}

// note: requires constructor T(string)
template <class T> 
void readListFile(list<T> &l, const char *name) {
  ifstream inputFile(name);
  if(!inputFile) 
    fatalError((string("Cannot open file ")+name));

  string aux;
  while( getline( inputFile , aux) )
    l.push_back(T(aux));

  inputFile.close();
}

template <class T> 
void readListFile(string sName, list<T> &l) {
  readListFile(l,sName.c_str());
}

template <class T> 
void readLinedFileWithBackInsert(T &l, const char *name) {
  back_insert_iterator< T > ii(l);

  ifstream inputFile(name);
  if(!inputFile) 
    fatalError((string("Cannot open file ")+name));

  string aux;
  while( getline( inputFile , aux) )
    *ii++ = aux;

  inputFile.close();
}

template <class T> 
void writeListFile(const char *name, list<T> &l) {
  ofstream outputFile(name);
  copy(l.begin(),l.end(),ostream_iterator<T>(outputFile,"\n"));
  outputFile.close();
}

void printStateVector(const StateVector &);

void printCompactStateVector(const StateVector &);

void printStateVectorSet(set<StateVector> &);

/* returns a random integer in the range [lower,upper] (ie. includes upper) */
int getRandNum(int lower,int upper);

template <class _InputIter>
inline bool verifyNoDuplicates(_InputIter iIter, _InputIter last)
{
  bool bResult = true;
  _InputIter jIter;

  while(iIter!=last && bResult) {
    jIter = iIter;
    jIter++;
    while(jIter!=last && bResult) {
      bResult = (*iIter == *jIter);
      jIter++;
    }

    iIter++;
  }

  return bResult;
}


template<class T>
bool verifyNoDuplicates(Vector< T > &list)
{
  int i,j;
  bool bResult = true;
  for(i=1; i<list.size() && bResult; i++)
    for(int j=0; j<i && bResult; j++)
      if(list[i] == list[j])
	bResult = false;
  return bResult;
}

void readExpressionInputFile(const char *, bool, bool ,MatrixGroup & ,Vector< discrete > &);

void readExpressionInputFile(string ,bool, bool,MatrixGroup & ,Vector< discrete > &);

char *itoa(int);

// A more general form I should move to soon.
template<class T>
void readExpressionInputFile(const char *inputFileName,
			     bool bWraparound,
			     GMatrixGroup<T> &vMatrices)
{
  int i, j;

  Vector< Vector< T > > auxIn;
  readFile(auxIn, inputFileName);

  Matrix<T> M(auxIn[0].size(),
	      auxIn.size()+(bWraparound?1:0) ); //(nGenes X #pts+0/1)

  FOR(i,auxIn[0].size())
    FOR(j,auxIn.size()) 
      M[i][j] = auxIn[j][i]; 

  if(bWraparound) 
    FOR(i,auxIn[0].size())
      M[i][auxIn.size()] = auxIn[auxIn.size()-1][i];

  vMatrices.push_back(M);
}

template<class T>
void readExpressionInputFile(string sInputFileName,
			     bool bWraparound,
			     GMatrixGroup<T> &vMatrices)
{
  readExpressionInputFile(sInputFileName.c_str(), bWraparound, vMatrices);
}

void readConstraintsFile(const char *, Vector< map<int,int> >& );

void readConstraintsFile(string &, Vector< map<int,int> >& );

template <class T>
void printObject(ostream &s, const T &t) 
{
  s << t;
}

namespace std {

template <class T, class W>
ostream& operator<<(ostream &s, const pair<const T,W> &p)
{ 
  return s << p.first << "," << p.second << " "; 
}

}

template <class T>
std::ostream& operator<<(std::ostream &s, const map<T,T> &p)
{ 
  copy(p.begin(),p.end(),
       ostream_iterator< pair<const T,T> >(s," ") );
  return s;
}

template <class T, class W>
std::ostream& operator<<(std::ostream &s, const map<T,W> &p)
{ 
  copy(p.begin(),p.end(),
       ostream_iterator< pair<const T,W> >(s," ") );
  return s;
}

template <class T>
std::ostream& operator<<(std::ostream &s, const vector<T> &p)
{ 
  copy(p.begin(),p.end(),
       ostream_iterator< const T >(s," ") );
  return s;
}

template <class T>
std::istream& operator>>(std::istream &s, vector<T> &p)
{ 
  T aux;
  while(s) {
    s >> aux;
    p.push_back(aux);
  }
  return s;
}

template <class T>
std::ostream& operator<<(std::ostream &s, const list<T> &p)
{ 
  copy(p.begin(),p.end(),
       ostream_iterator< const T >(s," ") );
  return s;
}

template <class T>
std::ostream& operator<<(std::ostream &s, const set<T> &p)
{ 
  copy(p.begin(),p.end(),
       ostream_iterator< const T >(s," ") );
  return s;
}

/*
template <class T>
void readLineFile(T<string> &l, const char *name) {
  ifstream inputFile(name);
  if(!inputFile) 
    fatalError("cannot open file ");

  string aux;
  while( getline( inputFile , aux) )
    l.push_back(aux);

  inputFile.close();
}
*/

template <class T>
void normalise( GMatrixGroup<T> &G )
{
  int i,j,k;
  T min,max;

  assert(G.size() && G[0].size());

  FOR(j,G[0].num_rows()) {
    min = max = G[0][j][0];
    FOR(i,G.size()) {
      FOR(k,G[i].num_cols()) {
	if(G[i][j][k] < min)
	  min = G[i][j][k];
	else if(G[i][j][k] > max)
	  max = G[i][j][k];
      }
    } // at this point min and max have been computed
    
    FOR(i,G.size()) {
      FOR(k,G[i].num_cols()) 
	G[i][j][k] = (G[i][j][k]-min)/(max-min);
    }

  }
}

template <class T>
void normalise( vector<T> &v )
{
  int i;
  T min,max;

  assert(v.size());
  min = max = v[0];
  FOR(i,v.size()) {
    if(v[i] < min)
      min = v[i];
    else if(v[i] > max)
      max = v[i];
  }

  if(max!=min) {
    FOR(i,v.size()) 
      v[i] = (v[i]-min)/(max-min);
  }    
}

/*
template <class T>
bool operator==(const vector<T> &A, 
		const vector<T> &B)
{
  bool bResult = true;

  Subscript i, M = A.size();
  if(M!=B.size())
    bResult = false;

  for (i=0; i<M && bResult; i++)
    bResult = (A[i] == B[i]);

  return bResult;
}

template <class T>
bool operator<(const vector<T> &A, 
	       const vector<T> &B)
{
  return lexicographical_compare(A.begin(),A.end(),B.begin(),B.end());
}
*/

template <class T>
void printTimeSeries( int x, vector< Matrix<T> > &M )
{
  int i,j;
  FOR(i,M.size()) {
    FOR(j,M[i].num_cols())
      cout << M[i][x][j] << " ";
    cout << "\t";
  }
  cout << endl;
}

template <class T>
Vector<T> getColumn(Matrix<T> M, int col)
{
  int i;
  assert((col>=0) && (col<M.num_cols()));
  Vector<T> vResult = Vector<T>(M.num_rows());
  FOR(i,M.num_rows())
    vResult[i] = M[i][col];
  return vResult;
}

template <class T>
Vector<T> getLastColumn(Matrix<T> M)
{
  return getColumn(M,M.num_cols()-1);
}

template <class T>
Vector<T> getRow(Matrix<T> M, int r)
{
  int i;
  assert((r>=0) && (r<M.num_rows()));
  Vector<T> vResult = Vector<T>(M.num_cols());
  FOR(i,M.num_cols())
    vResult[i] = M[r][i];
  return vResult;
}

// is a trajectory (in Matrix representation) in steady state or not?
bool isSteadyState(Matrix<int> &);

// returns number of links in a network
int parsimony(Matrix<int> &);

/*
template <class T>
bool operator==(const std::pair<T,T> &A,
		const std::pair<T,T> &B)
{
  return ( (A.first==B.first) && (A.second==B.second) );
}
*/

// TODO: template this ...
bool operator==(const std::map<int,int> &A,
		const std::map<int,int> &B);
/*
template <class T>
int sign(T x)
{
  int r = 0;
  if(x) {
    if(x<0)
      r = -1;
    else 
      r = 1;
  }
  cerr << "r=" << r << endl;
  return r;
}
*/
template <class X, class Y>
struct FirstLess : binary_function<pair<X,Y>,pair<X,Y>,bool> {
  bool operator() (const pair<X,Y>& p, const pair<X,Y>& q) const
  {
    return p.first < q.first;
  }
};

template <class X, class Y>
struct FirstGreater : binary_function<pair<X,Y>,pair<X,Y>,bool> {
  bool operator() (const pair<X,Y>& p, const pair<X,Y>& q) const
  {
    return p.first > q.first;
  }
};

template <class X, class Y>
struct SecondLess : binary_function<pair<X,Y>,pair<X,Y>,bool> {
  bool operator() (const pair<X,Y>& p, const pair<X,Y>& q) const
  {
    return p.second < q.second;
  }
};

template <class X, class Y>
struct SecondGreater : binary_function<pair<X,Y>,pair<X,Y>,bool> {
  bool operator() (const pair<X,Y>& p, const pair<X,Y>& q) const
  {
    return p.second > q.second;
  }
};

// these 3 from popsim_utils
bool fileIsDir(string sName);
void createCleanDirectory(string sDirName);
void createDirectory(string sDirName);

#define PR(x) cerr << #x << "="<< x << endl;

#define TB(x) x << "\t"

#define IFOR(it,obj) for(it=obj.begin();it!=obj.end();it++)

#define SQ(x) (x)*(x)

template <class U>
double 
sumSquareDifferences(const U &v1, const U &v2)
{
  double difference, sumSquareDifferences = 0.0;
  int i,N = v1.size();
  FOR(i,N) {
    difference = v1[i]-v2[i];
    sumSquareDifferences += (difference*difference);
  }
  sumSquareDifferences /= N;
  return sumSquareDifferences;
}

vector<double>
linspace(double lower, double upper, double stepSize);

Vector<double>
convertToTntVector(const gsl_vector *v);

#endif





