#include "util.h"
#include "BivariateNormalConditional.h"
#include "LineReader.h"
#include <gsl/gsl_cdf.h>

BivariateNormalConditional::
BivariateNormalConditional()
{
}

BivariateNormalConditional::
~BivariateNormalConditional()
{
}

void
BivariateNormalConditional::
getConditionalsV2(double x2, double &condMean, double &condSd)
{
  double mu1 = X.mean();
  double mu2 = Y.mean();
  double s1  = X.stddev();
  double s2  = Y.stddev();
  double rho = corrcoeff();
  condMean = mu1 + (s1/s2)*rho*(x2-mu2);
  condSd = sqrt( (1-pow(rho,2.0))*pow(s1,2.0) );
}

double
BivariateNormalConditional::
getV1ConditionalOnV2(double x1, double x2)
{
  double condMean, condSd, r;
  getConditionalsV2(x2,condMean,condSd);
  r = gsl_cdf_gaussian_P(x1-condMean,condSd);
  //PR( gsl_cdf_gaussian_P(x1-condMean,condSd) );
  //exit(0);
  return(r);
}

int BIVmain(int argc, char *argv[]) {
  BivariateNormalConditional bnc;
  ifstream testInputFile("/home/tom/tmp/bivarTest.dat");
  LineReader<double> lineReader(testInputFile);
  vector<double> vBuffer;
  while(lineReader.next(vBuffer)) {
    //cout << TB(vBuffer[0]) << vBuffer[1] << endl;
    bnc.sample(vBuffer[0],vBuffer[1]);
  }
  double condMean, condSd;
  bnc.getConditionalsV2(2.0,condMean,condSd);
  PR(condMean);
  PR(condSd);
  PR( bnc.getV1ConditionalOnV2(-3.0,2.0) );
}
