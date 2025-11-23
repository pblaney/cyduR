#include "MotifIdentifier.h"

MotifIdentifier::
MotifIdentifier()
{
}

MotifIdentifier::
MotifIdentifier(const Motif &_motif)
  :motif(_motif)
{
}

MotifIdentifier::
~MotifIdentifier()
{
}

bool
MotifIdentifier::
check(int m)
{
  const set<int> &sAux = motifReference[motif[state]];
  bool bResult = ( sAux.find(m) != sAux.end() ) ; // true if it's there
  return bResult;
}

// NOTE: inefficient?
vector<int>
MotifIdentifier::
find(const Sequence &sequence)
{
  int i;
  state = 0;
  vFound.clear();
  vFound.reserve(sequence.size()-motif.size());
  assert(sequence.size()>=motif.size());
  int lastPos = (sequence.size()-motif.size())+1;
  FOR(i,lastPos) { 
    state = 0;
    while( (state<motif.size()) && check(sequence[i+state]) ) {
      state++;
    }
    if(state == motif.size()) {
      vFound.push_back(i);
    }
  }
  return vFound;
}

vector<int>
MotifIdentifier::
shiftedFind(const Sequence &_sequence)
{
  int i;
  vector<int> vFound = find(_sequence);
  FOR(i,vFound.size())
    vFound[i] += motif.mutationPosition;
  return vFound;
}

