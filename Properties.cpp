#include "Properties.h"
#include "util.h"

#include <iostream>
#include <algorithm>

/**
 * Inserts default values into multimap. 
 * Empty virtual implementation called at start of load(). 
 * Can be overridden to insert defaults in subclasses.
 **/
void Properties::insertDefaults() {
}

/** 
 * Processes a line of text from input file.
 * @param sAux A line of text.
 **/
void Properties::processLine(string sAux) {
  int iPos1, iPos2;

  if( (sAux.size()>0) && (sAux[0]=='#') )
    return;

  iPos1 = sAux.find_first_of("=");

  if( iPos1 < sAux.size() ) {    
    string sKey = sAux.substr(0,iPos1);
    iPos1++;
    
    do {
	iPos2 = sAux.find_first_of(":",iPos1);
	string sValue = sAux.substr(iPos1,iPos2-iPos1);
	insert(make_pair(sKey,sValue));
	iPos1 = iPos2+1;
    } while( iPos2 < sAux.size() );

  }
}

/**
 * Load from a stream into the multimap.
 * @param isIn An input stream to load from
 **/
void Properties::load(istream &isIn) {
  insertDefaults();

  string sAux;
  while( getline( isIn , sAux) )
    processLine(sAux);
}

/** 
 * Get a single property from the multimap.
 * Retrieves a single value from the multimap for the key.
 * There may be several values, and the function will retrieve only one.
 * If no value is found, the default is returned
 * @param sKey key to look for. 
 * @param sDefault default value to return if no value found.
 **/ 
string Properties::getProperty(string sKey, string sDefault) {
  iterator itAux = find(sKey);
  if( itAux != end() ) 
    return itAux->second;
  else
    return sDefault;
}

/** 
 * Get a single property from the multimap.
 * Retrieves a single value from the multimap for the key.
 * There may be several values, and the function will retrieve only one.
 * If no value is found, an empty string ("") is returned
 * @param sKey key to look for. 
 **/ 
string Properties::getProperty(string sKey) {
  return getProperty(sKey,"");
}

/**
 * Get a single property from the multimap.
 * See function above.
 **/
string Properties::getProperty(string sKey, const char *cDef) {
  return getProperty(sKey,string(cDef));
}

/** 
 * Get list of properties from the multimap.
 * If no value is found, an empty vector is returned
 * @param sKey key to look for. 
 **/ 
vector<string> Properties::getPropertyList(string sKey) {
  vector<string> vResult;
  Properties::iterator itAux = find(sKey);
  while( (itAux!=end()) && (itAux->first==sKey))
    vResult.push_back( (itAux++)->second );
  return vResult;
}

/**
 * Boolean check for a specific key/value pair.
 * @param sKey Key to check for.
 * @param sValue Value to check for
 * @return true if key/value pair exists, false otherwise.
 **/

bool Properties::checkKeyValue(string sKey, string sValue) {
  bool bResult = false;
  typedef const_iterator I;
  pair<const_iterator,const_iterator> b = equal_range(sKey);
  for(const_iterator i=b.first; i!=b.second && !bResult; ++i) 
    bResult = (i->second == sValue);
  return bResult;
}

/**
 * Boolean check for a specific key/value pair.
 * Same as function above
 **/
bool Properties::checkKeyValue(const char *cKey, const char *cValue) {
  return checkKeyValue(string(cKey),string(cValue));
}


