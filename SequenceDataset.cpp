#include "SequenceDataset.h"

SequenceDataset::
SequenceDataset()
{
}

SequenceDataset::
~SequenceDataset()
{
}

SequenceDataset operator+(const SequenceDataset &a,
			  const SequenceDataset &b)
{
  SequenceDataset c(a);
  c.vvSeq = concatenate(a.vvSeq,b.vvSeq);
  return c;
}

