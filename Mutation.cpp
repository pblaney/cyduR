#include "Mutation.h"
#include "util.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
extern gsl_rng *r_gsl_rng;

Mutation::
Mutation()
  : isNull(false)
{
}

Mutation::
~Mutation()
{
}

void
Mutation::
add(const set<int> &_from, const set<int> &_to)
{
  from.push_back(_from);
  to.push_back(_to);
}

bool
Mutation::
check(int baseA, int baseB)
{
  int i;
  bool bResult = false;
  for(i=0; (i<from.size()) && (!bResult); i++) {
    if( ( from[i].find(baseA) != from[i].end()) && ( to[i].find(baseB) != to[i].end() ) ) 
      bResult = true;
  }
  return bResult;
}

bool
Mutation::
apply(Sequence &vAux, int pos, double pChange)
{
  int i;
  Mutation &auxMutation = *this;
}
