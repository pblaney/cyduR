#include "util.h"
#include "DiscreteHistogram.h"

DiscreteHistogram::DiscreteHistogram()
  :dh(1,(state)0) {
}

DiscreteHistogram::~DiscreteHistogram() {
}

void DiscreteHistogram::sample(const discrete &x) {
  int i;
  if(x >= dh.size()) { //need to expand
    StateVector oldDh(dh);
    dh.newsize(x+1);
    for(i=0; i<oldDh.size(); i++)
      dh[i] = oldDh[i];
    for(i=oldDh.size(); i<dh.size(); i++)
      dh[i] = 0;
  }
  dh[x]++;
}

discrete DiscreteHistogram::getLowerBound() {
  return 0;
}

// Note: returns -1 if no samples presented
discrete DiscreteHistogram::getMinSampled() {
  int i, pos=-1;
  bool found = false;

  for(i=0; i<dh.size() && !found; i++)
    if(dh[i]) {
      found = true;
      pos = i;
    }

  return pos;
}

discrete DiscreteHistogram::getUpperBound() {
  return dh.size();
}

discrete DiscreteHistogram::operator[] (const discrete &i) {
  if(i<dh.size())
    return dh[i];
  else
    return 0;
}

discrete DiscreteHistogram::getSum() {
  int i,sum=0;
  FOR(i,dh.size())
    sum += dh[i];
  return sum;
}
