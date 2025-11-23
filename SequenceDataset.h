#ifndef SequenceDataset_h
#define SequenceDataset_h

#include "util.h"
#include "Sequence.h"

class SequenceDataset 
{
public:
  SequenceDataset();
  ~SequenceDataset();
  string label;
  Sequence vRef;
  Vector< Sequence > vvSeq;
};

SequenceDataset operator+(const SequenceDataset &a,
			  const SequenceDataset &b);

#endif
