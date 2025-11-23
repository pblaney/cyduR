#ifndef ContinuousHistogram_h
#define ContinuousHistogram_h

#include <cmath>

#include "global.h"
#include "DiscreteHistogram.h"

class ContinuousHistogram {
protected:
  double min, max, stepsize;
  int nBuckets;
  DiscreteHistogram dh;
public:
  ContinuousHistogram(double min, double max, int nBuckets);
  ~ContinuousHistogram();

  template <class T> void sample(T _x) {
    double x = _x;
    if(x>=min && x<=max) {
      int dumpBucket = (int)floor((x-min)/stepsize);
      dh.sample(dumpBucket);
    }
    else
      cerr << "ContinuousHistogram error: sample outside range ignored \t" 
	   << x << endl;
  };
  
  discrete getNBuckets();
  double getStepsize();

  discrete getLowerBound();
  discrete getUpperBound();
  discrete getMinSampled();
  discrete operator[](const discrete &);
}; 

#endif
