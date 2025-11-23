#ifndef BivariateNormalConditional_h
#define BivariateNormalConditional_h

#include "StatsSampler2Vars.h"

class BivariateNormalConditional : public StatsSampler2Vars
{
public:
  BivariateNormalConditional();
  ~BivariateNormalConditional();
  void getConditionalsV2(double x2, double &condMean, double &condSd); // get mean of variable 1 conditional on variable 2 = x2
  double getV1ConditionalOnV2(double x1, double x2);
};

#endif
