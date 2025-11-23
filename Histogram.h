#ifndef Histogram_h
#define Histogram_h

// virtual Histogram class

#include "global.h"
#include "StateVector.h"

class Histogram {
 public:
  virtual void sample(const discrete &) =0;
  virtual discrete getLowerBound() =0;
  virtual discrete getUpperBound() =0;
  virtual discrete operator[](const discrete &) =0;
};

#endif
