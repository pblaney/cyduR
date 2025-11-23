#include "util.h"
#include <cstdlib>
#include <cstdio>
#include <iostream>

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <unistd.h>
#include <errno.h>

extern int errno;


// TODO: implement base conversions for ultra-efficient StateVector representations

void fatalError(char *errMessage) {
  cerr << errMessage << endl;
  exit(-1);
}

void fatalError(string errMessage) {
  cerr << errMessage << endl;
  exit(-1);
}

void printStateVector(const StateVector &s) {
    for(int i=0; i<s.size(); i++) {
      cout.width(4);
      cout << s[i];
    }
    cout << endl;
}

void printCompactStateVector(const StateVector &s)
{
  for(int i=0; i<s.size(); i++) {
    if( s[i] < 0 )
      cout << "-";
    else if( s[i] > 0 )
      cout << "+";
    else 
      cout << "0";
  }
  cout << endl;
}

void printStateVectorSet(set< StateVector > &targetSet) {
  int i;
  set< StateVector >::iterator it = targetSet.begin();
  while(it != targetSet.end()) {
    const StateVector &s = *it++;
    for(i=0; i<s.size(); i++) {
      cout.width(4);
      cout << s[i];
    }
    cout << endl;
  }
}

/* returns a random integer in the range [lower,upper] (ie. includes upper) */
int getRandNum(int lower, int upper) {
    double range = (upper - lower) + 1;
    return ( lower + (int) (range*(rand()/(RAND_MAX+1.0))) );
}

/*
void readListFile(list<string> &l, const char *name) {
  ifstream inputFile(name);
  if(!inputFile) 
    fatalError("cannot open file ");

  string aux;
  while( getline( inputFile , aux) )
    l.push_back(aux);

  inputFile.close();
}
*/

void readExpressionInputFile(const char *inputFileName
			     ,bool definesDefaults
			     ,bool bWraparound
			     ,MatrixGroup &vMatrices
			     ,Vector< discrete > &vDefaults) 
{
  int i, j;

  Vector< Vector< discrete > > auxIn;
  readFile(auxIn, inputFileName);
  if(definesDefaults)
    vDefaults = auxIn[0];

  Matrix<discrete> M(auxIn[0].size(),
		     auxIn.size()+(bWraparound?1:0) ); //(nGenes X #pts+0/1)

  FOR(i,auxIn[0].size())
    FOR(j,auxIn.size()) 
      M[i][j] = auxIn[j][i]; 

  if(bWraparound) 
    FOR(i,auxIn[0].size())
      M[i][auxIn.size()] = auxIn[auxIn.size()-1][i];

  vMatrices.push_back(M);
}

void readExpressionInputFile(string sInputFileName
			     ,bool definesDefaults
			     ,bool bWraparound
			     ,MatrixGroup &vMatrices
			     ,Vector< discrete > &vDefaults) 
{
  readExpressionInputFile(sInputFileName.c_str()
			  ,definesDefaults,bWraparound,vMatrices,vDefaults);
}

char *itoa(int x)
{
  static char buf[20];
  sprintf(buf,"%d",x);
  return buf;
}

void readConstraintsFile(const char *inputFileName, 
			 Vector< map<int,int> > &vConstraints)
{
  int i, j;
  int POSITION=0, VALUE=1; // two columns within each matrix

  Vector< Matrix<int> > vmConstraints;
  readFile(vmConstraints,inputFileName); 

  vConstraints.newsize(vmConstraints.size());
  
  FOR(i,vConstraints.size())
    FOR(j,vmConstraints[i].num_rows())
    vConstraints[i].insert(make_pair(vmConstraints[i][j][POSITION],
				     vmConstraints[i][j][VALUE]));
  /*  
  map<int,int>::iterator cur;
  FOR(i,vConstraints.size()) {
    for(cur=vConstraints[i].begin(); cur!=vConstraints[i].end(); cur++)
      cerr << (*cur);
    cerr << endl;
  }
  */
}

void readConstraintsFile(string &sInputFileName, 
			 Vector< map<int,int> > &vConstraints)
{
  readConstraintsFile(sInputFileName.c_str(), vConstraints);
}

bool isSteadyState(Matrix<int> &mResult) {
  assert(mResult.num_cols()>=2);
  int numStates = mResult.num_cols();
  if(!(getColumn(mResult,numStates-1)==getColumn(mResult,numStates-2)))
      return false;
  return true;
}

int parsimony(Matrix<int> &C) {
  int i,j,nGenes=C.num_rows(),result=0;
  assert(nGenes=C.num_cols());
  FOR(i,nGenes)
    FOR(j,nGenes)
      if(C[i][j])
	result++;
  return result;
}

bool operator==(const std::map<int,int> &A,
		const std::map<int,int> &B)
{
  int i;
  bool bResult = true;
  if(A.size() != B.size())
    return false;
  
  map<int,int>::const_iterator at = A.begin();
  map<int,int>::const_iterator bt = B.begin();
  while(at!=A.end() && bResult) {
    bResult = (at->first==bt->first) && (at->second==bt->second);
    at++;
    bt++;
  }

  return bResult;
}


bool fileIsDir(string sName)
{
  bool bResult = false;
  DIR *dp = opendir(sName.c_str());
  if(dp != NULL) {
    bResult = true;  
    closedir(dp);
  }
  return bResult;
}

void createCleanDirectory(string sDirName)
{
  DIR *existingDirectory = opendir(sDirName.c_str());
  if(existingDirectory == NULL) {
    //cerr << "creating directory: " << sDirName << endl;
    if(errno == ENOENT) // directory does not exist ...
      mkdir(sDirName.c_str(),0755);
    else
      assert(false);
  }
  else {
    struct dirent *dirp;
    while( (dirp=readdir(existingDirectory)) != NULL ) {
      string sDName = dirp->d_name;
      string sEntryName = sDirName + "/" + sDName;
      if(sDName != "." && sDName!="..") {
	if(!fileIsDir(sEntryName)) {
	  //cerr << "deleting file:" << sEntryName << endl;
	  unlink(sEntryName.c_str());
	}
	else {
	  //cerr << "deleting subdir:" << sEntryName << endl;
	  createCleanDirectory(sEntryName); // USE WITH CARE!
	  rmdir(sEntryName.c_str());
	}
      }
    }
    closedir(existingDirectory);
  }
}

void createDirectory(string sDirName)
{
  DIR *existingDirectory = opendir(sDirName.c_str());
  if(existingDirectory == NULL) {
    //cerr << "creating directory: " << sDirName << endl;
    if(errno == ENOENT) // directory does not exist ...
      mkdir(sDirName.c_str(),0755);
    else
      assert(false);
  }
  else
    closedir(existingDirectory);
}

vector<double>
linspace(double lower, double upper, double stepSize)
{
  int i;
  vector<double> vL;
  for(i=0; i==0 || (*vL.rbegin()<upper); i++) {
    vL.push_back(lower+stepSize*i);
  }
  return vL;
}

Vector<double>
convertToTntVector(const gsl_vector *v)
{
  int i;
  Vector<double> r( v->size );
  FOR(i,v->size)
    r[i] = gsl_vector_get(v,i);
  return r;
}
