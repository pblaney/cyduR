#include "RandomizedIota.h"

#include <numeric>
#include <algorithm>
using namespace std;
/*
template <class _ForwardIterator, class _Tp>
void 
iota(_ForwardIterator __first, _ForwardIterator __last, _Tp __value)
{
  while (__first != __last)
    *__first++ = __value++;
}
*/
RandomizedIota::
RandomizedIota(int n)
  :vector<int>(n)
{
  iota(begin(),end(),0);
  shuffle();
}

RandomizedIota::
~RandomizedIota()
{
}

void
RandomizedIota::
shuffle()
{
  random_shuffle(begin(),end());
}

