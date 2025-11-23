#include "MotifMutationFrequency.h"

MotifMutationFrequency::
MotifMutationFrequency()
{
}

MotifMutationFrequency::
MotifMutationFrequency(Motif _motif, Mutation _mutation, double _freq) 
{ 
  motif=_motif; 
  mutation=_mutation; 
  freq=_freq;
}

MotifMutationFrequency::
~MotifMutationFrequency()
{
}  

