#ifndef Properties_h
#define Properties_h

#include <string>
#include <map>
#include <set>
#include <vector>

using namespace std;

/**
 * A utility class for configuration files.
 * Reads and writes configuration files in a format similar to 
 * the Java Properties class. The text file format is:
 * key=value1:value2:etc.
 * Inherits STL multimap of string to string
 */

class Properties : public multimap<string,string> {
 protected:
  virtual void insertDefaults();
  void processLine(string);
 public:
  void load(istream &);
  string getProperty(string, string);
  string getProperty(string);
  vector<string> getPropertyList(string sKey);
  string getProperty(string, const char *);
  bool checkKeyValue(string,string);
  bool checkKeyValue(const char *,const char *);
};


#endif
