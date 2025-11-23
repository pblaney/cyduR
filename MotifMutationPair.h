#ifndef MotifMutationPair_h
#define MotifMutationPair_h

#include "Motif.h"
#include "Mutation.h"

class MotifMutationPair
{
public:
  MotifMutationPair();
  MotifMutationPair(Motif _motif, Mutation _mutation);
  ~MotifMutationPair();
  Motif motif;
  Mutation mutation;
};

#endif
