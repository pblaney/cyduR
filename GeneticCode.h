#ifndef GeneticCode_h
#define GeneticCode_h

#include "util.h"

class GeneticCode
{
public:
  GeneticCode();
  ~GeneticCode();

  map<string,char> codonMap;
  map<char,vector<string> > aaMap;

  void addEntry(char *_codon, char *_aa);
};

#endif
