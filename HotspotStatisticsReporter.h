#ifndef HotspotStatisticsReporter_h
#define HotspotStatisticsReporter_h

#include "util.h"
#include "Sequence.h"
#include "GeneticCode.h"
#include "MotifReference.h"
#include "MotifIdentifier.h"
#include "FrequencyDependentRandomizer.h"

struct CombinedMotif 
{
public:
  int bSenseStrandOnly;
  MotifIdentifier fwdMI, revMI;
  string label;
  string getPrettyString();
};

class HotspotStatisticsReporter
{
public:
  enum {WRC,SYC,CCC,TC,CG,TTC,TCC,YC,BC,S_TA,S_TT};
  enum {BOTH_STRANDS,SENSE_STRAND_ONLY};

  GeneticCode geneticCode;
  MotifReference motifReference;

  MotifIdentifier wrcMI, gywMI, sycMI, grsMI, cccMI, gggMI, tcMI, gaMI, cgMI, ttcMI, gaaMI, tccMI, ggaMI, ycMI, grMI, bcMI, gvMI ,s_taMI, s_ttMI;

  vector<CombinedMotif> vCM;

  map< char, FrequencyDependentRandomizer* > mapFdr;

  bool bAssignedNucleotideProportions;
  Vector<double> vProportionsACG;

  HotspotStatisticsReporter();
  ~HotspotStatisticsReporter();

  void addSearchMotif(string sMotif, int mutatedPosition, int bSenseStrand);
  void assignDefaultMotifs();
  void readExternalMotifs(ifstream &_motifStream);

  double calcGcProportion(Sequence sOriginal);
  void computeCodonWeights(double gcProportion);

  void evaluate( map<string,string> &mss, bool bCoding);

  char getAA(const Sequence &sAux, int ntPosition);
  pair<char,string> getAAwithCodon(const Sequence &sAux, int ntPosition);

  Sequence randomizeCodons(const Sequence &sAux);
  Vector<int> countHotspots(Sequence &sAux);

  void assignNucleotideComposition(double aProportion, double cProportion, double gProportion);

  void countReplacements(Sequence &sAux, MotifIdentifier &mi, int mutationTo, int &numReplacements, int &numHotspots);
  Matrix<int> countTransitionReplacements(Sequence &sAux);
  
  string changeNANtoNA(double x);
};

#endif
