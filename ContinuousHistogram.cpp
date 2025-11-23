#include "ContinuousHistogram.h"

ContinuousHistogram::
ContinuousHistogram(double _min, double _max, int _nBuckets) 
  : min(_min), max(_max), nBuckets(_nBuckets)
{
  assert(max>=min);
  stepsize = (max-min)/nBuckets;
}  

ContinuousHistogram::
~ContinuousHistogram()
{
}

discrete ContinuousHistogram::getNBuckets()
{
  return nBuckets;
}

double ContinuousHistogram::getStepsize()
{
  return stepsize;
}

discrete ContinuousHistogram::
getLowerBound()
{
  return dh.getLowerBound();
}

discrete ContinuousHistogram::
getUpperBound()
{
  return dh.getUpperBound();
}

discrete ContinuousHistogram::
getMinSampled()
{
  return dh.getMinSampled();
}

discrete ContinuousHistogram::operator[](const discrete &i)
{
  return dh[i];
}
