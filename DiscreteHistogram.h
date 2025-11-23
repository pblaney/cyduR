#ifndef DiscreteHistogram_h
#define DiscreteHistogram_h

#include "global.h"
#include "Histogram.h"
#include "StateVector.h"

class DiscreteHistogram :public Histogram {
 protected:
  discrete baseOffset;
  WeightVector dh;
 public:
  DiscreteHistogram();
  ~DiscreteHistogram();
  void sample(const discrete &);
  discrete getLowerBound();
  discrete getMinSampled();
  discrete getUpperBound();
  discrete getSum();
  discrete operator[](const discrete &);
};

#endif
