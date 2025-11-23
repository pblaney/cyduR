#ifndef RandomizedIota_h
#define RandomizedIota_h

#include <vector>

using namespace std;

class RandomizedIota: public vector<int>
{
public:
  RandomizedIota(int);
  ~RandomizedIota();  
  void shuffle();
};

#endif
