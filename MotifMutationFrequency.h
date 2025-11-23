#ifndef MotifMutationFrequency_h
#define MotifMutationFrequency_h

#include "Motif.h"
#include "Mutation.h"

class MotifMutationFrequency
{
public:
  MotifMutationFrequency();
  MotifMutationFrequency(Motif _motif, Mutation _mutation, double _freq);
  ~MotifMutationFrequency();
  Motif motif;
  Mutation mutation;
  double freq;
};

#endif
