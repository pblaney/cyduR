#ifndef Motif_h
#define Motif_h

#include "tnt.h"

class Motif : public Vector<int>
{
public:
  int mutationPosition;
  Motif();
  Motif(int size, int defaultValue, int mutationPosition);
  Motif(int size, int mutationPosition);
  ~Motif();
  string getString() const;
  string getPrettyString() const; 
};

#endif
