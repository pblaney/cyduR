#ifndef LineReader_h
#define LineReader_h

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <strstream>

template <class T>
class LineReader {
 protected:
  istream &inputStream;
  string buffer;
  vector<T> tBuffer;

  void parseBuffer() { // parse buffer into stateBuffer
    istrstream ist(buffer.c_str());
    T aux;
    tBuffer.clear();
    while(ist >> aux) {
      tBuffer.push_back(aux);
    }
  }

 public:
  LineReader(istream &_inputStream)
    :inputStream(_inputStream)
  {
    if(!inputStream) 
      fatalError("cannot open data input stream");
  }
  
  ~LineReader() {}

  bool next(vector<T> &param){
    if(getline(inputStream,buffer)) {
      parseBuffer();
      param = tBuffer;
    }
    if(inputStream.eof())
      return false;
    else
      return true;
  }
};

#endif



