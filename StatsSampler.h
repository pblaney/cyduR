#ifndef StatsSampler_h
#define StatsSampler_h

#include <cmath>
#include <functional>

class StatsSampler {
protected:
  double min, max;
  double sum, sum_squares;
  unsigned long long int n;

public:
  StatsSampler();
  ~StatsSampler();
  void initialise();

  template <class T> void sample(T x) {
    double _x = (double)x;
    sum += _x;
    sum_squares += ( _x * _x );

    if(min == HUGE_VAL)
      min = max = _x;
    else if(_x < min)
      min = _x;
    else if(_x > max)
      max = _x;

    n++;
  };

  double getMin();
  double getMax();
  double getSum();
  double getSumSquares();
  double getN();
  double mean() const; 
  double var();
  double stddev();
};

#endif
