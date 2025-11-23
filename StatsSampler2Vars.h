#ifndef StatsSampler2Vars_h
#define StatsSampler2Vars_h

#include "StatsSampler.h"

class StatsSampler2Vars {
protected:
  double sum_products;
public:
  StatsSampler X, Y;
  StatsSampler2Vars();
  ~StatsSampler2Vars();
  void initialise();

  template <class U, class V> void sample(U x, V y) {
    double _x = (double)x;
    double _y = (double)y;
    X.sample( _x );
    Y.sample( _y );
    sum_products += _x*_y;
  };

  double cov();
  double corrcoeff();
};

#endif
