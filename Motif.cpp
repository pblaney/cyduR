#include "util.h"
#include "Motif.h"

Motif::
Motif()
{
}

Motif::
Motif(int _size, int _default, int _mutationPosition)
  : Vector<int>(_size,_default), mutationPosition(_mutationPosition)
{
}

Motif::
Motif(int _size, int _mutationPosition)
  : Vector<int>(_size,14), mutationPosition(_mutationPosition)
{
}

Motif::
~Motif()
{
}

string
Motif::
getString() const 
{  
  const Motif &motif = *this;
  int j;
  const char *mutations[] = {"G","C","A","T","S","R","K","M","Y","W","H","D","B","V","N"};
  string r;
  FOR(j,size()) {
    r += mutations[motif[j]];
  }
  return r;
}

string
Motif::
getPrettyString() const
{  
  const Motif &motif = *this;
  int j;
  const char *mutations[] = {"G","C","A","T","S","R","K","M","Y","W","H","D","B","V","N"};
  string r;
  bool bPending = false;
  FOR(j,size()) {
    if(bPending) {
      r += "_";
      bPending = false;
    }
    if(j==mutationPosition) {
      r += "_";
      bPending = true;
    }
    r += mutations[motif[j]];
  }
  if(bPending)
    r += "_";
  return r;
}
