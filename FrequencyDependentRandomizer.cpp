#include <cmath>

#include "util.h"
#include "FrequencyDependentRandomizer.h"

gsl_rng *FrequencyDependentRandomizer::r = NULL;

void
FrequencyDependentRandomizer::
setSeed(int seed)
{
  if(r==NULL) {
    r = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(r,seed);
  }
}

FrequencyDependentRandomizer::
FrequencyDependentRandomizer()
  :grdt(NULL)
{
  setSeed(0);
}

FrequencyDependentRandomizer::
FrequencyDependentRandomizer(const Vector<double> &_vF)
  :grdt(NULL)
{
  if(r==NULL) 
    r = gsl_rng_alloc(gsl_rng_default);
  setFrequencies(_vF);
}


FrequencyDependentRandomizer::
~FrequencyDependentRandomizer()
{
  if(grdt!=NULL) gsl_ran_discrete_free(grdt);
  //gsl_rng_free(r); 
}

void
FrequencyDependentRandomizer::
setFrequencies(const Vector<double> &vF)
{
  int i, N=vF.size() ;
  double *ptrP = new double[N];
  FOR(i,N)
    ptrP[i] = vF[i];

  if(grdt!=NULL) gsl_ran_discrete_free(grdt);
 
  grdt = gsl_ran_discrete_preproc(N,ptrP);

  delete [] ptrP;
}

int 
FrequencyDependentRandomizer::
getRand()
{
  int result = gsl_ran_discrete(r,grdt);
  return result;
}

// the number of independent trials that need to be performed 
// so that the probability of missing a category (category prob.=1/n)
// is exactly k (defined as a static constant in header).
/*
int 
FrequencyDependentRandomizer::
calculateK(int n)
{
  if(n<=1) return 1;
  double result = log(z)/(log((double)n-1)-log((double)n));
  result = ceil(result);
  return (int)result;
}
*/

double expfun(int _k, double K)
{
  //double K = 2.222;
  double k = (double) _k;
  double r = (1.0-exp(-1.0/K)) * exp(-k/K);
  return r;
}

double PLfun(int _k)
{
  double k = (double) _k;
  double K = 2.295;
  double r = ( (1/k) * (exp(-k/K)) ) / (-log(1.0-exp(-1.0/K)));
  return r;
}

#include "StatsSampler.h"

double processExp()
{
  int i,j;
  int max_k = 50;
  double desiredExponent = 0.45;
  double desiredK = 1.0/desiredExponent;
  
  FrequencyDependentRandomizer fdr1;
  Vector<double> f1(max_k,0.0);

  // note need for difference between expfun and PLfun (which needs shift +1)
  FOR(i,max_k)
    f1[i] = expfun(i,desiredK);

  fdr1.setFrequencies(f1);

  StatsSampler ss;
  FOR(j,1000000) {
    int x = fdr1.getRand();
    ss.sample(x);
  }

  return(ss.mean());
}

double processPL()
{
  int i,j;
  int max_k = 50;
  
  FrequencyDependentRandomizer fdr1;
  Vector<double> f1(max_k,0.0);

  // note need for difference between expfun and PLfun (which needs shift +1)
  FOR(i,max_k)
    f1[i] = PLfun(i+1);

  fdr1.setFrequencies(f1);

  StatsSampler ss;
  FOR(j,1000000) {
    int x = fdr1.getRand();
    ss.sample(x+1);
  }

  return(ss.mean());
}
/*
int main(int argc, char *argv[])
{
  cerr << "Mean using exponential ..." << processExp() << endl;
  cerr << "Mean using power law ..." << processPL() << endl;
}
*/






