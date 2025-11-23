#include <math.h>
#include <iostream>
using namespace std;

#include "StatsSampler2Vars.h"

StatsSampler2Vars::StatsSampler2Vars()
{
  initialise();
}

StatsSampler2Vars::~StatsSampler2Vars()
{
}

void StatsSampler2Vars::initialise()
{
  X.initialise();
  Y.initialise();
  sum_products = (double)0.0;
}

double StatsSampler2Vars::cov()
{
  return (sum_products/X.getN()) - (X.mean()*Y.mean());
}

double StatsSampler2Vars::corrcoeff()
{
  return cov() / (X.stddev()*Y.stddev());
}


/*
int main(int argn, char *argv[]) {
  int i;
  int X[] = {80,50,36,58,72,60,56,68};
  int Y[] = {65,60,35,39,48,44,48,61};

  StatsSampler ss1;
  StatsSampler2Vars ss2;

  for(i=0;i<8;i++)
    ss2.sample(X[i],Y[i]);
  cout << "cov=" << ss2.cov() << "\tcorrcoeff=" << ss2.corrcoeff() << "\n" ;
}
*/
