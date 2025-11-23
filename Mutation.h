#ifndef Mutation_h
#define Mutation_h

#include "util.h"
class Sequence;

class Mutation
{
public:
  Mutation();
  ~Mutation();
  bool isNull;
  vector< set<int> > from, to;
  void add(const set<int> &from, const set<int> &to);
  bool check(int baseA, int baseB);
  bool apply(Sequence &vAux, int pos, double pChange);
};

#endif
