#ifndef TYPEDEFS_H
#define TYPEDEFS_H

/*
  Global typedefs definitions
*/

#define ON true
#define OFF false

typedef int state;

typedef int discrete;

typedef double continuous;

enum message {CONTINUE, STUCK, ERROR, FINAL} ;

using namespace std;

#define FOR(k,nIter) for(k=0;k<nIter;k++) 
#define MAX(x,y) ((x>y)?x:y) 
#define MIN(x,y) ((x<y)?x:y) 

#endif
