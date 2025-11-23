#ifndef Sequence_h
#define Sequence_h

#include "tnt.h"
#include "util.h"
#include "Mutation.h"

class Sequence : public Vector<int>
{
public:
  Sequence();
  Sequence(int, int);
  Sequence( Vector<int> );
  ~Sequence();
  //string sLabel;
  int distance(const Sequence &, bool bIgnoreGaps=false) const;
  Matrix<int> mutationCounts(const Sequence &y) const;
  Vector<int> baseCounts() const;
  int distance(const Sequence &, const Mutation &);
  static int getComplement(int x);
  Sequence complement();
  Sequence cut(int l, int r) const;
  Sequence paste(const Sequence &y);
  void print() const;
  void print(ostream &outFile) const;
};

#endif
