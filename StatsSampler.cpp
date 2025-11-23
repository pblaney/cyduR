#include "StatsSampler.h"

StatsSampler::StatsSampler() 
{
  initialise();
}

StatsSampler::~StatsSampler() 
{
}

void StatsSampler::initialise()
{
  min = max = HUGE_VAL;
  sum = sum_squares = 0;
  n = 0;
}

double StatsSampler::getMin() {
  return min;
}

double StatsSampler::getMax() {
  return max;
}

double StatsSampler::getSum() {
  return sum;
}

double StatsSampler::getSumSquares() {
  return sum_squares;
}

double StatsSampler::getN() {
  return (double)n;
}

double StatsSampler::mean() const {
  return sum/(double)n;
}

double StatsSampler::var() {
  return (sum_squares - ((sum*sum)/(double)n )) / (double)n;
}

double StatsSampler::stddev() {
  return sqrt(var());
}
