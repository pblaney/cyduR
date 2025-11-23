#include "scharff_utils.h"
#include "generate.h"
#include "Sequence.h"
#include "MotifReference.h"

Sequence::
Sequence()
{
}

Sequence::
Sequence(int _length, int _x)
  :Vector<int>(_length, _x)
{
}

Sequence::
Sequence(Vector<int> _x)
  :Vector<int>(_x)
{
}

Sequence::
~Sequence()
{
}

int
Sequence::
distance(const Sequence &y, bool bIgnoreGaps) const
{
  int i;
  const Sequence &x = *this;
  int differenceCount = 0;
  if(bIgnoreGaps) {
    FOR(i,x.size())
      if( (x[i]!=y[i]) && (x[i]<4) && (y[i]<4) )
	differenceCount++;
  } else {
    FOR(i,x.size())
      if(x[i]!=y[i])
	differenceCount++;
  }
  return differenceCount;
}

Matrix<int>
Sequence::
mutationCounts(const Sequence &y) const
{
  int i;
  const Sequence &x = *this;
  Matrix<int> mResult(5,5,0);
  FOR(i,x.size())
    mResult[x[i]][y[i]]++;
  return mResult;
}

Vector<int>
Sequence::
baseCounts() const
{
  int i;
  const Sequence &x = *this;
  Vector<int> vCounts(5,0); // might be a DiscreteHistogram
  FOR(i,x.size())
    vCounts[ x[i] ]++;
  return vCounts;
  
}

// This is the distance taking into account a specific mutation
// Note that it is not commutative: this=mutation.from, operand=mutation.to
int
Sequence::
distance(const Sequence &y, const Mutation &mutation)
{
  int i,j;
  const Sequence &x = *this;
  int differenceCount = 0;
  FOR(i,x.size()) {
    FOR(j,mutation.from.size()) {
      if( (mutation.from[j].find(x[i])!=mutation.from[j].end()) && 
	  (mutation.to[j].find(y[i])!=mutation.to[j].end()) )
	differenceCount++;
    }
  }
  return differenceCount;
}

int 
Sequence::
getComplement(int x)
{
  int r;
  switch(x) {
  case G: r=C; break;
  case C: r=G; break;
  case A: r=T; break;
  case T: r=A; break;
  case GAP: r=GAP; break;
  case U: r=A; break; // G or A?
  default: PR(x); assert(false); break;
  }
  
  return r;
}

Sequence
Sequence::
complement()
{
  int i;
  Sequence &seqThis = *this;
  Sequence result(seqThis);
  FOR(i,result.size())
    result[i] = getComplement(result[i]);
  return result;
}

Sequence
Sequence::
cut(int l, int r) const
{  
  const Sequence &x = *this;
  int i, j, length=MAX( 0, (r-l)+1 );
  Sequence result(length,0);
  for(i=l, j=0; i<=r; i++, j++)
    result[j] = x[i];
  return result;
}

Sequence
Sequence::
paste(const Sequence &y)
{
  Sequence &x = *this;
  Sequence result = concatenate_vector(x,y);
  return result;
}

void 
Sequence::
print() const
{
 int i; 
 const Sequence &seqThis = *this;
 FOR(i,seqThis.size())
   cout << MotifReference::bases[ seqThis[i] ];
}

void 
Sequence::
print(ostream &outFile) const
{
 int i; 
 const Sequence &seqThis = *this;
 FOR(i,seqThis.size())
   outFile << MotifReference::bases[ seqThis[i] ];
}

