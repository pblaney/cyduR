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
#include "generate.h"

typedef map<int,int> svType;
typedef set< pair<int,int> > setLinkType;

//#define DEBUG

gsl_rng *r_gsl_rng = NULL;

void rng_initialise(int seed)
{
  if(r_gsl_rng == NULL) {
    r_gsl_rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(r_gsl_rng,seed);
  }
}

void rng_free()
{
  if(r_gsl_rng != NULL)
    gsl_rng_free(r_gsl_rng);  
}

void
reportLogLogDegreeDist(const Matrix <discrete> &C)
{
  int i,j;
  DiscreteHistogram dh;

  PR(C);
  
  FOR(i,C.num_rows()) {
    int degree=0;
    FOR(j,C.num_cols())
      if( (C[i][j] != 0) || (C[j][i] != 0) )
	degree++;
    dh.sample(degree);
  }
  
  for(i=1; i<=dh.getUpperBound(); i++) {
    double logk = log((double)i);
    double logPk = log((double)dh[i]);
    if(logPk >= 0.0)
      cout << logk << "\t" << logPk << endl;
  }
}


vector<double>
getDegreeDist(const Matrix <discrete> &C)
{
  int i,j;
  DiscreteHistogram dh;
  
  FOR(i,C.num_rows()) {
    int degree=0;
    FOR(j,C.num_cols())
      if( C[i][j] || C[j][i] )
	degree++;
    dh.sample(degree);
  }

  int sum = 0;
  for(i=dh.getLowerBound(); i<=dh.getUpperBound(); i++)
    sum += dh[i];

  vector<double> r(dh.getLowerBound()+dh.getUpperBound()+1,0.0);
  FOR(i,r.size())
    r[i] = (double)dh[i]/sum;
  return r;
}

void generateP(int N, Vector<double> &P, Matrix<discrete> &C)
{
  int i,j, sum=0;
  N = P.size();

  FOR(i,N) {
    P[i] = 0.0;
    FOR(j,N) {
      if(C[i][j] || C[j][i]) { // count both inbound and outbound
	P[i] += 1.0;
	sum++;
      }
    }
  }
}

Matrix<discrete>
createBarabasiNetwork(int maxT)
{
  int i,j,k, m0=1;//, maxT = 10;
  Matrix<discrete> C(maxT,maxT,(discrete)0);
  rng_initialise(0);

  for(i=m0; i<maxT; i++) {
    int newSize = i+1;
    Vector<double> P(newSize);
    generateP(newSize,P,C);

    double *ptrP = new double[newSize];
    FOR(j,newSize)
      ptrP[j] = P[j];

    gsl_ran_discrete_t *grdt = gsl_ran_discrete_preproc(newSize,ptrP);
    int randomNode = gsl_ran_discrete(r_gsl_rng,grdt);
    C[i][randomNode] = C[randomNode][i] = 1; // undirected graph

    delete [] ptrP;
    gsl_ran_discrete_free(grdt);
  }
  return C;
}

Matrix<discrete>
createDirectedBarabasiNetwork(int maxT) {
  Matrix<discrete> C = createBarabasiNetwork(maxT);
  convertUndirectedNetwork(C);
  return C;
}

void
convertUndirectedNetwork(Matrix<discrete> &C)
{
  int i,j,N=C.num_rows();
  FOR(i,N) {
    FOR(j,i+1) {
      if(C[i][j] != 0) {
	int randomPick = getRandNum(0,1);
	int randomSign = getRandNum(0,1)*2-1;
	if(randomPick) {
	  C[i][j] = randomSign;
	  C[j][i] = 0;
	} else {
	  C[i][j] = 0;
	  C[j][i] = randomSign;	  
	}
      }
    }
  }
}

Matrix<discrete>
createErdosRenyiNetwork(int N, double p) {
  int i,j;
  Matrix<discrete> C(N,N,0);
  rng_initialise(0);
  FOR(i,N) {
    for(j=0; j<=i; j++) {
      if(gsl_ran_flat(r_gsl_rng,0.0,1.0) < p) {
	if(getRandNum(0,1))
	  C[i][j] = getRandNum(0,1)*2-1;
	else
	  C[j][i] = getRandNum(0,1)*2-1;
      }
    }
  }
  return C;
}

Matrix<discrete>
createUnsignedErdosRenyiNetwork(int N, double p) {
  int i,j;
  Matrix<discrete> C(N,N,0);
  rng_initialise(0);
  FOR(i,N) {
    for(j=0; j<=i; j++) {
      if(gsl_ran_flat(r_gsl_rng,0.0,1.0) < p) {
	if(getRandNum(0,1))
	  C[i][j] = 1;
	else
	  C[j][i] = 1;
      }
    }
  }
  return C;
}

double 
averageConnectivity(const Matrix<discrete> &C)
{
  int i,j, sum=0, N=C.num_rows();
  FOR(i,N) {
    FOR(j,N)
      if(C[i][j] != 0)
	sum++;
  }
  return (double)sum / (N*N);  
}
/*
int
countNonZeroEntries(const Matrix<discrete> &C)
{
  int i,j, sum=0, N=C.num_rows();
  FOR(i,N) {
    FOR(j,N)
      if(C[i][j] != 0)
	sum++;
  }
  return sum;
}

vector< pair<int,int> >
extractNonZeroEntries(const Matrix<discrete> &C)
{
  vector< pair<int,int> > vR;
  int i,j, sum=0, N=C.num_rows();
  FOR(i,N) {
    FOR(j,N)
      if(C[i][j] != 0)
	vR.push_back( make_pair(i,j) );
  }
  return vR;
}
*/
Matrix<discrete>
createUnsignedSiegalBergmanNetwork(int numGenes,double c)
{
  int i,j;
  Matrix<discrete> mGeneInteractions(numGenes,numGenes);

  FOR(i,numGenes) {
    FOR(j,numGenes) {
      if( gsl_ran_flat(r_gsl_rng,0.0,1.0) < c )
	mGeneInteractions[i][j] = 1;
    }
  }
  return mGeneInteractions;
}






