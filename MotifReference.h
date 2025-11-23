#ifndef MotifReference_h
#define MotifReference_h

#include "util.h"

class MotifReference : public vector< set<int> >
{
public:
  MotifReference();
  ~MotifReference();

  static char bases[];
  static int reverseCodes[];

  map<char,int> mBP;

  void add(char *_s, int pos);

  static bool isBase(char);

};

#endif
