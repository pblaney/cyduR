#include "util.h"
#include "Sequence.h"
#include "ContinuousHistogram.h"
#include <gsl/gsl_vector.h>
#include "EgnetProperties.h"
#include "SequenceDataset.h"

enum {G,C,A,T,GAP,U,EMPTY};
enum {SENSE,ANTISENSE};

double calcDifference(const vector<int> &a, const vector<int> &b);
double calcDifference(const Vector<int> &a, const Vector<int> &b);

Sequence
convertSequence(const string &sData, int maxPosition=-1);

Vector< Sequence > 
convertSequence(const map<string,string> &mss, int maxPosition=-1);

void 
writeFastaSequence(string sName, Sequence sequence, ostream &outFile);

void
writeFastaSequences(Sequence vRef, Vector<Sequence> vSeqs, string sFilename);

vector<int>
extractPolymorphicSites( Sequence vRef, Vector<Sequence> vvSeq );

int
countPolymorphicSites( Sequence vRef, Vector<Sequence> vvSeq );

void
removeGapsInReferenceSequence( Sequence &vRef, Vector<Sequence> &vvSeq );

void
removeUnmutatedSequences( const Sequence &vRef, Vector<Sequence> &vvSeq );

int
countMutations(const Sequence &vRef, Vector<Sequence> &vvSeq );

void
chopSequences( Sequence &vRef, Vector<Sequence> &vvSeq, int maxPosition );

void
chopSequences( Vector<Sequence> &vvSeq, int maxPosition );

void
chopSequence( Sequence &vRef, int maxPosition );

void
prechopSequences( Sequence &vRef, Vector<Sequence> &vvSeq, int cutoff );

Vector<double>
buildSiteMutationProfile(Sequence vRef, Vector<Sequence> vvSeq);

Mutation
createMutation(Vector<char> vC);

Mutation
createMutationToAny(char c);

ContinuousHistogram
buildProfile(Sequence vRef, Vector<Sequence> vvSeq);

double 
convertLambda(double lambda, Sequence seqConsensus);

string
removeLeadingPath(string str);

SequenceDataset
readSingleDataset(string sInputFile, string sRefKey);

vector<SequenceDataset>
readMultipleDatasets(EgnetProperties propsConfig, string p_sRefKey=string("refKey"), string p_sInputFiles=string("inputFiles") );

Matrix<int> 
generateDotMatrix(const SequenceDataset &sequenceDataset);

double
extractMeanMutability(const Vector<double> &vMutationProfile,
		      const vector<int> &vPositions);

template <class T>
Matrix<double> 
normalize(const Matrix<T> &x)
{
  int i,j,sum;
  Matrix<double> r(x.num_rows(),x.num_cols(),0.0);
  FOR(i,x.num_rows()) {
    FOR(j,x.num_cols()) 
      sum += x[i][j];
  }

  if(sum) {
    FOR(i,x.num_rows()) {
      FOR(j,x.num_cols())
	r[i][j] = (double)x[i][j] / sum;
    }
  }
  
  return r;
}

template <class T>
T sgn(T x)
{
  T r = 0;
  if(x<0)
    r = -1;
  else {
    if(x>0)
      r = 1;
  }
  return r;
}

template <class T>
Matrix<T> sgn(const Matrix<T> &A)
{
  int i,j;
  Matrix<int> B(A);

  for(i=0;i<B.num_rows();i++) {
    for(j=0;j<B.num_cols();j++)
      B[i][j] = sgn(B[i][j]);
  }
  return B;
}

template <class T>
vector<T> sgn(const vector<T> &A)
{
  int i,j;
  vector<T> B(A);
  for(i=0;i<B.size();i++) 
    B[i] = sgn(B[i]);
  return B;
}

void 
multipleAssign(Vector<double> &x, vector<int> y, double z);

enum {HOT,NEUTRAL,COLD,OVERALL};

Vector<double>
calculateMutability(SequenceDataset sd,
		    const Vector<double> &vMutationProfile, 
		    const Vector<int> &vMutationsPerSequence,
		    bool bTopStrand);

Vector<double>
calculateMutabilityWithSiteCounts(SequenceDataset sd,
				  const Vector<double> &vMutationProfile, 
				  const Vector<int> &vMutationsPerSequence,
				  bool bTopStrand,
				  Vector<int> &vSiteCounts,
				  Vector<int> &vSites);

/*
template <class U, class V>
void multipleAssign(TNT::Vector<U> &x)//, vector<int> y, U z)
{
  int i;
  FOR(i,y.size())
    x[ y[i] ] = z;
}
*/
