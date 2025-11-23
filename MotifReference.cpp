#include <ctype.h>
#include "scharff_utils.h"
#include "MotifReference.h"

char MotifReference::bases[] = {'G','C','A','T','-'};
int MotifReference::reverseCodes[] = {1,0,3,2,4,8,7,6,5,9,11,10,13,12,14};

MotifReference::
MotifReference()
{
  int i;
  vector< set<int> > &vThis = *this;
  const char *IUPAC[] ={"G","C","A","T", "S", "R", "K", "M", "Y", "W",  "H",  "D",  "B",  "V",   "N"};
  const char *sRep[] = {"G","C","A","T","GC","GA","GT","CA","CT","AT","CAT","GAT","GCT","GCA","GCAT"};
        //               0   1   2   3    4    5    6    7    8    9    10    11    12    13     14

  //                     C   G   T   A   CG   CT   CA   GT   GA   TA   GTA   CTA   CGA   CGT   CGTA   reverse letters
  //                     1   0   3   2    4    8    7    6    5    9    11    10    13    12     14
        



  mBP['G'] = G;
  mBP['C'] = C;
  mBP['A'] = A;
  mBP['T'] = T;
  for(i=4; i<15; i++)   // adding IUPAC codes: 30/OCT/2014
    mBP[ IUPAC[i][0] ] = i;
  
  vThis.resize(15);
  FOR(i,15)
    add((char *)sRep[i],i);  
}

MotifReference::
~MotifReference()
{
}

void
MotifReference::
add(char *_s, int pos)
{
  int i;
  string s(_s);
  set<int> sAux;
  FOR(i,s.size())
    sAux.insert(mBP[s[i]]);

  vector< set<int> > &vThis = *this;
  vThis[pos] = sAux;
}

bool 
MotifReference::
isBase(char x)
{
  bool result = false;
  int ux = toupper(x);
  if( ux=='G' || ux=='C' || ux=='A' || ux=='T' )
    result = true;
  return result;
}
