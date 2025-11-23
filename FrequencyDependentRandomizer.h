#ifndef FrequencyDependentRandomizer_h
#define FrequencyDependentRandomizer_h

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "tnt.h"

class FrequencyDependentRandomizer
{
protected:
  static gsl_rng *r;
  gsl_ran_discrete_t *grdt;
public:
  //static double z = 0.25;
  static void setSeed(int seed);
  FrequencyDependentRandomizer();
  FrequencyDependentRandomizer(const Vector<double> &vF);
  ~FrequencyDependentRandomizer();

  void setFrequencies(const Vector<double> &vF);  
  int getRand();
  
  //static int calculateK(int);
};

#endif
