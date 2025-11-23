#ifndef SimpleFastaReader_h
#define SimpleFastaReader_h

class SimpleFastaReader
{
public:
  map<string,string> read(const char *fName);
  
  SimpleFastaReader();
  ~SimpleFastaReader();
};

#endif
