#ifndef MotifIdentifier_h
#define MotifIdentifier_h

#include "util.h"
#include "Motif.h"
#include "Sequence.h"
#include "MotifReference.h"

class MotifIdentifier
{
protected:
  int state;
  vector<int> vFound;
public:
  Motif motif;
  MotifIdentifier();
  MotifIdentifier(const Motif &);
  ~MotifIdentifier();
  MotifReference motifReference;

  vector<int> find(const Sequence &sequence);
  vector<int> shiftedFind(const Sequence &sequence);
  inline bool check(int);
};

#endif
