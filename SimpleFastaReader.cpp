#include "util.h"
//#include "stress_utils.h"

#define IFOR(it,obj) for(it=obj.begin();it!=obj.end();it++)

#include "SimpleFastaReader.h"
#include "LineReader.h"
#include "DiscreteHistogram.h"
#include "EgnetProperties.h"

SimpleFastaReader::
SimpleFastaReader()
{
}

SimpleFastaReader::
~SimpleFastaReader()
{
}

map<string,string>
SimpleFastaReader::
read(const char *fName)
{
  int i;
  map<string,string> mapResult;

  vector<string> vBuffer;
  ifstream inFile(fName);
  if(!inFile) {
    cerr << "could not open file " << fName << endl;
    exit(EXIT_FAILURE);
  }
  LineReader<string> lr0(inFile);
  string sKey, sData;
  bool bFirstTime = true;
  while(lr0.next(vBuffer)) {
    if( vBuffer.size() && vBuffer[0].size() ) { // i.e. line not empty
      if(vBuffer[0][0] == '>') {
	if(!bFirstTime) { // need to write results into the map
	  mapResult[sKey] = sData;
	  sData = "";
	}
	sKey = vBuffer[0];
	for(i=1; i<vBuffer.size(); i++)
	  sKey += " "+vBuffer[i];
	bFirstTime = false;
      }
      else
	sData += vBuffer[0];
    }
  }
  if(sKey.size())
    mapResult[sKey] = sData; // write last line in ...
  inFile.close();
  return mapResult;
}

int SFRmain(int argc, char *argv[])
{
  int i,j;

  if(argc < 2)
    fatalError("usage: egnet <config_file>");
  ifstream configFile( argv[1] );

  if(!configFile)
    fatalError("cannot open configuration file ...");
  EgnetProperties propsConfig;
  propsConfig.load( configFile );
  configFile.close();

  EgnetProperties::iterator itAux;

  //int maxPosition = 420;
  int maxPosition = atoi(propsConfig.getProperty("maxPosition").c_str());
  int comparisonKey = atoi(propsConfig.getProperty("comparisonKey").c_str());
  string sInputFile = propsConfig.getProperty("inputFile");
  string sOutputDir = propsConfig.getProperty("outputDir");

  createDirectory(sOutputDir);

  SimpleFastaReader sfr;
  map<string,string> mss = sfr.read(sInputFile.c_str());

  vector<string> vRefKeys, vRefData;
  itAux = propsConfig.find("refKeys");
  while( (itAux!=propsConfig.end()) && (itAux->first=="refKeys"))
    vRefKeys.push_back( (itAux++)->second );

  FOR(i,vRefKeys.size()) {
    vRefData.push_back( mss[vRefKeys[i]] ); 
    mss.erase(vRefKeys[i]);
    assert(vRefData[i].size()>maxPosition);
  }
  int numSeqs = mss.size();

  enum {G,C,A,T,GAP};
  Matrix<int> mSeq(numSeqs,maxPosition,GAP);
  map<string,string>::iterator it;
  i=0;
  IFOR(it,mss){
    string sData = it->second;
    for(j=0; j<sData.size() && j<maxPosition; j++) {
      switch(sData[j]) {
      case 'G': mSeq[i][j] = G; break;
      case 'C': mSeq[i][j] = C; break;
      case 'A': mSeq[i][j] = A; break;
      case 'T': mSeq[i][j] = T; break;
      default: break; // keep GAP
      }
    }
    i++;
  }

  map<char,int> mBP;
  mBP['G'] = G;
  mBP['C'] = C;
  mBP['A'] = A;
  mBP['T'] = T;
  mBP['-'] = GAP;
  char mutations[] = {'G','C','A','T'};
  int mutation;

  FOR(mutation,4) {
    char mutationFrom = mutations[mutation];
    string sFileName("from.X.dat");
    sFileName[5] = mutationFrom;
    ofstream outFile((sOutputDir+sFileName).c_str());
    FOR(i,maxPosition) {
      DiscreteHistogram dh;
      if( vRefData[comparisonKey][i]==mutationFrom ) {
	FOR(j,numSeqs) {
	  if(mSeq[j][i] != mBP[vRefData[comparisonKey][i]]) // there has been a mutation
	    dh.sample(mSeq[j][i]);
	}
      }
      outFile << TB(i);
      FOR(j,4)
	outFile << TB((double)dh[j]/numSeqs);
      outFile << endl;
    }
    outFile.close();
  }
}

