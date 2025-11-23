#ifndef StateVector_h
#define StateVector_h

#include "tnt.h"
#include "global.h"
#include <set>
#include <string>

class StateVector: public Vector<state> {
 protected:
  bool same(const Vector<state>&, const Vector<state>&) const;
 public:
  StateVector();
  StateVector(Subscript, const state&);
  StateVector(Subscript, char *);
  StateVector(Subscript, const state *);
  StateVector(const char *);
  StateVector(Vector<state>&);
  bool operator==(const Vector<state>&) const;
  bool operator==(const StateVector&) const;
  bool operator!=(const StateVector&);
  bool operator<(const StateVector&) const;
  int distance(const StateVector&) const;
  string getSignedRepr();
  StateVector &operator++(int);
  StateVector &operator--(int);
  StateVector &operator-();
  typedef state* iterator;
  iterator begin();
  iterator end();
  StateVector& resize(Subscript N);

  template<class T>
  StateVector &operator=(const vector<T> &X) {
    StateVector &a = *this;
    a.resize(X.size());
    for(int i=0;i<a.size();i++) 
      a[i] = (discrete) X[i];
    return *this;
  }

};

typedef StateVector WeightVector;

#endif
