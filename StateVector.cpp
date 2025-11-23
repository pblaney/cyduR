#include "StateVector.h"

#include <cstdlib>
#include <iostream>

StateVector::StateVector() :Vector<state>() {}

StateVector::StateVector(Subscript N, const state& v) :Vector<state>(N,v) {}

StateVector::StateVector(Subscript N, char *s) :Vector<state>(N,s) {}

StateVector::StateVector(Subscript N, const state *s) :Vector<state>(N,s) {}

StateVector::StateVector(Vector<state> &v) :Vector<state>(v) {}

StateVector::StateVector(const char *cRep) :Vector<state>() 
{
  string sRep(cRep);
  newsize(sRep.length());
  Vector<state> &_this = *this;

  for(int i=0; i<size(); i++) {
    switch( sRep[i] ) {
    case '-':  _this[i] =-1; break;
    case '0':  _this[i] = 0; break;
    default :  _this[i] = 1; break;
    }
  }
}

string StateVector::getSignedRepr() {
  string sResult = "";
  Vector<state> &_this = *this;

  for(int i=0; i<size(); i++) {
    switch( _this[i] ) {
    case -1 :  sResult += '-' ; break;
    case  0 :  sResult += '0' ; break;
    default :  sResult += '+' ; break;
    }
  }
  return sResult;
}

inline bool StateVector::same(const Vector<state> &a, const Vector<state> &b) const {
  bool bSame = true;
  int i=0;

  if(a.size() != b.size())
    bSame = false;

  while((i<a.size()) && (bSame==true) ) {
    bSame = (a[i]==b[i]);
    i++;
  }
  return bSame;
}

bool StateVector::operator==(const Vector<state> &b) const
{ return same(*this,b); }

bool StateVector::operator==(const StateVector &b) const
{ return same(*this,b); }

bool StateVector::operator!=(const StateVector &b)
{ return !same(*this,b); }

bool StateVector::operator<(const StateVector &b) const {
  const StateVector &a = *this;
  int i=0;
  while((i<a.size()) && (a[i]==b[i]) ) 
    i++;
  if(i<a.size() && a[i]<b[i])
    return true;
  else
    return false;
}

int StateVector::distance(const StateVector &b) const {
  const StateVector &a = *this;
  int i=0, sum=0;
  while(i<a.size()) {
    sum += abs(a[i] - b[i]); 
    i++;
  }
  return sum;
}

StateVector &StateVector::operator++(int) {
  const StateVector &a = *this;
  state *curPtr = (state *)&a[0], *endPtr = (curPtr+a.size());
  while( curPtr != endPtr )
    (*curPtr++)++;

  return *this;
}

StateVector &StateVector::operator--(int) {
  const StateVector &a = *this;
  state *curPtr = (state *)&a[0], *endPtr = (curPtr+a.size());
  while( curPtr != endPtr )
    (*curPtr++)--;

  return *this;
}

StateVector &StateVector::operator-() {
  const StateVector &a = *this;
  state *curPtr = (state *)&a[0], *endPtr = (curPtr+a.size());
  while( curPtr != endPtr ) {
    *curPtr = -(*curPtr);
    curPtr++;
  }
  return *this;
}

StateVector::iterator StateVector::begin() {
  StateVector &a = *this;
  return &a[0];
}

StateVector::iterator StateVector::end() {
  StateVector &a = *this;
  return (&a[0])+a.size();
}

StateVector& StateVector::resize(Subscript N) {
  return (StateVector &)newsize(N);
}
