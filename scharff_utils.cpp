#include "scharff_utils.h"
#include "MotifReference.h"
#include "DiscreteHistogram.h"
#include "analysis.h"
#include "MotifIdentifier.h"
#include "SequenceDataset.h"
#include "SimpleFastaReader.h"
#include "StatsSampler.h"

#include <algorithm>

double calcDifference(const vector<int> &a, const vector<int> &b)
{
  int i, diffCount=0;
  FOR(i,a.size())
    if(a[i] != b[i])
      diffCount++;
  return ( (double)diffCount / a.size() );
}

double calcDifference(const Vector<int> &a, const Vector<int> &b)
{
  int i, diffCount=0;
  FOR(i,a.size())
    if(a[i] != b[i])
      diffCount++;
  return ( (double)diffCount / a.size() );
}

Sequence
convertSequence(const string &sData, int maxPosition)
{
  int j;

  if(maxPosition==-1)
    maxPosition = sData.size(); // no maxPosition given, assume 

  Sequence vSeq(maxPosition,GAP);

  for(j=0; j<sData.size() && j<maxPosition; j++) {
    switch(sData[j]) {
      case 'G': 
      case 'g': 
	vSeq[j] = G; 
	break;
      case 'C': 
      case 'c': 
	vSeq[j] = C; 
	break;
      case 'A': 
      case 'a': 	
	vSeq[j] = A; 
	break;
      case 'T': 
      case 't': 
	vSeq[j] = T; 
	break;
      default: 
	break; // keep GAP
    }
  }

  return vSeq;
}

Vector< Sequence > 
convertSequence(const map<string,string> &mss, int maxPosition)
{
  Vector< Sequence > vvSeq(mss.size());
  map<string,string>::const_iterator it;
  int i=0;
  IFOR(it,mss){
    string sData = it->second;
    vvSeq[i] = convertSequence(sData,maxPosition);
    i++;
  }
  return vvSeq;
}

void 
writeFastaSequence(string sName, Sequence sequence, ostream &outFile)
{
  int i;
  outFile << ">" << sName << endl;
  FOR(i,sequence.size())
    outFile << MotifReference::bases[sequence[i]];
  outFile << endl;
}

void
writeFastaSequences(Sequence vRef, Vector<Sequence> vSeqs, string sFilename)
{
  int i;
  ofstream outFile(sFilename.c_str());
  writeFastaSequence(string("CONSENSUS"),vRef,outFile);
  FOR(i,vSeqs.size()) {
    string sName = string("SEQ")+itoa(i);
    writeFastaSequence(sName,vSeqs[i],outFile);
  }
  outFile.close();
}

vector<int>
extractPolymorphicSites( Sequence vRef, Vector<Sequence> vvSeq )
{
  int i, j, numSites = vRef.size();
  vector<int> vResult;
  FOR(i,numSites) {
    bool bDifferenceFound = false;
    j=0;
    while(!bDifferenceFound && j<vvSeq.size()) {
      if(vvSeq[j][i] != vRef[i]) {
	vResult.push_back(i);
	bDifferenceFound = true;
      }
      j++;
    }
  }
  return vResult;
}

int
countPolymorphicSites( Sequence vRef, Vector<Sequence> vvSeq )
{
  return extractPolymorphicSites(vRef,vvSeq).size();
}

void
removeGapsInReferenceSequence( Sequence &vRef, Vector<Sequence> &vvSeq )
{
  int i,j;
  // 1. build a list of gaps in vRef
  vector<int> vGaps,vIndex;
  FOR(i,vRef.size()) {
    if(vRef[i] == GAP)
      vGaps.push_back(i);
    else 
      vIndex.push_back(i);
  }

  // 2. create the new ref sequence
  Sequence aux_vRef(vIndex.size(),0);
  FOR(i,aux_vRef.size())
    aux_vRef[i] = vRef[ vIndex[i] ];
  
  // 3. now cut the gaps out of the mutant sequences
  Vector<Sequence> aux_vvSeq(vvSeq.size(),Sequence(vIndex.size(),0));
  FOR(i,aux_vvSeq.size()) {
    FOR(j,aux_vvSeq[i].size()) 
      aux_vvSeq[i][j] = vvSeq[i][ vIndex[j] ];
  }

  vRef = aux_vRef;
  vvSeq = aux_vvSeq;
}

void
removeUnmutatedSequences( const Sequence &vRef, Vector<Sequence> &vvSeq )
{
  int i;
  list<int> mutatedIndices;
  FOR(i,vvSeq.size()) {
    if( vRef.distance(vvSeq[i]) > 0 )
      mutatedIndices.push_back(i);
  }
  list<int>::iterator it;
  Vector<Sequence> aux_vvSeq(mutatedIndices.size());
  i=0;
  IFOR(it,mutatedIndices)
    aux_vvSeq[i++] = vvSeq[*it];
  vvSeq = aux_vvSeq;
}

int
countMutations(const Sequence &vRef, Vector<Sequence> &vvSeq )
{
  int i, result = 0;
  FOR(i,vvSeq.size()) 
    result += vRef.distance(vvSeq[i]);
  return result;
}

void
chopSequences( Sequence &vRef, Vector<Sequence> &vvSeq, int maxPosition )
{
  int i,j;
  if(maxPosition < vRef.size()) {
    Sequence aux_vRef(maxPosition,0);
    FOR(i,maxPosition)
      aux_vRef[i] = vRef[i];

    Vector<Sequence> aux_vvSeq(vvSeq.size(),Sequence(maxPosition,0));
    FOR(i,aux_vvSeq.size()) {
      FOR(j,aux_vvSeq[i].size()) 
	aux_vvSeq[i][j] = vvSeq[i][j];
    }
    
    vRef = aux_vRef;
    vvSeq = aux_vvSeq;
  }
}

void
chopSequences( Vector<Sequence> &vvSeq, int maxPosition )
{
  int i,j;
  if( (maxPosition>=0) && (vvSeq.size()>0) && (maxPosition<vvSeq[0].size()) ) {
    Vector<Sequence> aux_vvSeq(vvSeq.size(),Sequence(maxPosition,0));
    FOR(i,aux_vvSeq.size()) {
      FOR(j,aux_vvSeq[i].size()) 
	aux_vvSeq[i][j] = vvSeq[i][j];
    }
    vvSeq = aux_vvSeq;
  }
}

void
chopSequence( Sequence &vRef, int maxPosition )
{
  int i,j;
  if(maxPosition < vRef.size()) {
    Sequence aux_vRef(maxPosition,0);
    FOR(i,maxPosition)
      aux_vRef[i] = vRef[i];
    vRef = aux_vRef;
  }
}

void
prechopSequences( Sequence &vRef, Vector<Sequence> &vvSeq, int cutoff )
{
  int i,j;
  int maxPosition = vRef.size()-cutoff;

  if(maxPosition < vRef.size()) {
    Sequence aux_vRef(maxPosition,0);
    FOR(i,maxPosition)
      aux_vRef[i] = vRef[i+cutoff];

    Vector<Sequence> aux_vvSeq(vvSeq.size(),Sequence(maxPosition,0));
    FOR(i,aux_vvSeq.size()) {
      FOR(j,maxPosition) 
	aux_vvSeq[i][j] = vvSeq[i][j+cutoff];
    }
    
    vRef = aux_vRef;
    vvSeq = aux_vvSeq;
  }
}

Vector<double>
buildSiteMutationProfile(Sequence vRef, Vector<Sequence> vvSeq)
{
  int i,numSequences = vvSeq.size();
  Matrix<DiscreteHistogram> mDh = generateInternalMutationHistograms(vRef,vvSeq);
  Vector<double> vResult(vRef.size());

  FOR(i,vRef.size()) {
    int siteSum = mDh[vRef[i]][i].getSum();
    vResult[i] = ((double)siteSum/numSequences);
  }
  return vResult;
}

Mutation
createMutation(Vector<char> vC)
{
  MotifReference motifReference;
  set<int> from, to;
  from.insert(motifReference.mBP[vC[0]]);
  to.insert(motifReference.mBP[vC[1]]);
  Mutation mutation;
  mutation.add(from,to);
  return mutation;
}

Mutation
createMutationToAny(char c)
{
  int i;
  MotifReference motifReference;
  set<int> from, to;
  from.insert(motifReference.mBP[c]);
  FOR(i,4) {
    if(c != motifReference.bases[i])
      to.insert(motifReference.mBP[ motifReference.bases[i] ]);
  }
  Mutation mutation;
  mutation.add(from,to);
  return mutation;
}

ContinuousHistogram
buildProfile(Sequence vRef, Vector<Sequence> vvSeq)
{
  int i,numSequences = vvSeq.size();

  ContinuousHistogram ch(0,1,50);
  Matrix<DiscreteHistogram> mDh = generateInternalMutationHistograms(vRef,vvSeq);

  FOR(i,vRef.size()) {
    int siteSum = mDh[vRef[i]][i].getSum();
    if(siteSum) // there was at least one mutation
      ch.sample((double)siteSum/numSequences);
  }
  return ch;
}

double 
convertLambda(double lambda, Sequence seqConsensus)
{
  Sequence seqComplement = seqConsensus.complement();
  int L = seqConsensus.size()*2;
  int numWrc=0, numC=0;

  {
    Motif motif(3,2); motif[0] = 9; motif[1] = 5; motif[2] = 1; // WRC
    MotifIdentifier motifIdentifier(motif); 
    numWrc += motifIdentifier.shiftedFind(seqConsensus).size();
  }

  {
    Motif motif(3,2); motif[2] = 9; motif[1] = 5; motif[0] = 1; // CRW
    MotifIdentifier motifIdentifier(motif); 
    numWrc += motifIdentifier.shiftedFind(seqComplement).size();
  }

  {
    Motif motif(1,0); motif[0] = 1; // C
    MotifIdentifier motifIdentifier(motif); 
    numC += motifIdentifier.shiftedFind(seqConsensus).size();
    numC += motifIdentifier.shiftedFind(seqComplement).size();
  }

  double pMotif = (lambda*numWrc) / ( (L-numWrc) + (lambda*numWrc) );
  double pNonMotifC = (1-pMotif) * ( (numC-numWrc) / (double)(L-numWrc) );
  double pRatio = pMotif / pNonMotifC;
  /*
  PR(numWrc);
  PR(L);
  PR(lambda);
  PR(pMotif)
  PR(numC)
  PR(pNonMotifC)
  PR(pRatio);
  */
  return pRatio;
}

string
removeLeadingPath(string str)
{
  int pos = str.find_last_of("/");
  if(pos != string::npos)
    str = str.substr(pos+1);
  return str;
}

SequenceDataset
readSingleDataset(string sInputFile, string sRefKey)
{
  SequenceDataset sdResult;
  SimpleFastaReader sfr;
  map<string,string> sData = sfr.read(sInputFile.c_str());
  string sRef = sData[sRefKey];
  sData.erase(sRefKey);
  sdResult.vRef = convertSequence(sRef);
  sdResult.vvSeq = convertSequence(sData);
  sdResult.label = removeLeadingPath(sInputFile);
  return sdResult;
}

vector<SequenceDataset>
readMultipleDatasets(EgnetProperties propsConfig, string p_sRefKey, string p_sInputFiles)
{
  vector<SequenceDataset> vResult;
  string sRefKey = propsConfig.getProperty(p_sRefKey);

  EgnetProperties::iterator itAux = propsConfig.find(p_sInputFiles);
  while( (itAux!=propsConfig.end()) && (itAux->first==p_sInputFiles)) {
    string sInputFile = (itAux++)->second;
    SequenceDataset sdAux = readSingleDataset(sInputFile,sRefKey);
    vResult.push_back(sdAux);
  }
  return vResult;
}

Matrix<int> 
generateDotMatrix(const SequenceDataset &sequenceDataset)
{
  int i,j;

  Matrix<int> mResult(sequenceDataset.vvSeq.size(),sequenceDataset.vRef.size(),0);

  FOR(i,sequenceDataset.vvSeq.size()) {
    FOR(j,sequenceDataset.vRef.size()) {
      int ref = sequenceDataset.vRef[j];
      int seq = sequenceDataset.vvSeq[i][j];
      if( (ref != seq) && (ref<4) && (seq<4) )
	mResult[i][j] = 1;
    }
  }
  return mResult;
}

double
extractMeanMutability(const Vector<double> &vMutationProfile,
		      const vector<int> &vPositions)
{
  int i;
  StatsSampler ss;
  FOR(i,vPositions.size())
    ss.sample(vMutationProfile[ vPositions[i] ]);
  return(ss.mean());
}


void 
multipleAssign(Vector<double> &x, vector<int> y, double z)
{
  int i;
  FOR(i,y.size())
    x[ y[i] ] = z;
}

Vector<double>
calculateMutability(SequenceDataset sd,
		    const Vector<double> &vMutationProfile, 
		    const Vector<int> &vMutationsPerSequence,
		    bool bTopStrand)
{
  Vector<Sequence> vResult;
  int i;

  Motif wrcMotif(3,2); wrcMotif[0] = 9; wrcMotif[1] = 5; wrcMotif[2] = C; // WRC
  Motif sycMotif(3,2); sycMotif[0] = 4; sycMotif[1] = 8; sycMotif[2] = C; // SYC
  Motif cMotif(1,0);  cMotif[0] = C; // just C

  Motif gywMotif(3,0); gywMotif[0] = G; gywMotif[1] = 8; gywMotif[2] = 9; // GYW
  Motif grsMotif(3,0); grsMotif[0] = G; grsMotif[1] = 5; grsMotif[2] = 4; // GRS
  Motif gMotif(1,0);  gMotif[0] = G; // just G

  Motif hotMotif, coldMotif, allMotif;
  if(bTopStrand) {
    hotMotif = wrcMotif;
    coldMotif = sycMotif;
    allMotif = cMotif;
  } else {
    hotMotif = gywMotif;
    coldMotif = grsMotif;
    allMotif = gMotif;
  }

  MotifIdentifier allMotifIdentifier(allMotif);
  vector<int> allPositions = allMotifIdentifier.shiftedFind(sd.vRef);

  MotifIdentifier hotMotifIdentifier(hotMotif);
  vector<int> hotPositions = hotMotifIdentifier.shiftedFind(sd.vRef);

  MotifIdentifier coldMotifIdentifier(coldMotif);
  vector<int> coldPositions = coldMotifIdentifier.shiftedFind(sd.vRef);

  set<int> setAux, setNeutrals;
  set_difference(allPositions.begin(), allPositions.end(), 
		 hotPositions.begin(), hotPositions.end(),
                 inserter(setAux, setAux.begin()) );

  set_difference(setAux.begin(), setAux.end(), 
		 coldPositions.begin(), coldPositions.end(),
                 inserter(setNeutrals, setNeutrals.begin()) );
  
  vector<int> neutralPositions(allPositions.size()-(hotPositions.size()+coldPositions.size()));
  copy(setNeutrals.begin(),setNeutrals.end(),neutralPositions.begin());

  double hotMutability = extractMeanMutability(vMutationProfile,hotPositions);
  double coldMutability = extractMeanMutability(vMutationProfile,coldPositions);
  double neutralMutability = extractMeanMutability(vMutationProfile,neutralPositions);
  double overallMutability = extractMeanMutability(vMutationProfile,allPositions);

  Vector<double> vResults(4);
  vResults[HOT] = hotMutability;
  vResults[COLD] = coldMutability;
  vResults[NEUTRAL] = neutralMutability;
  vResults[OVERALL] = overallMutability;

  return(vResults);
}

Vector<double>
calculateMutabilityWithSiteCounts(SequenceDataset sd,
				  const Vector<double> &vMutationProfile, 
				  const Vector<int> &vMutationsPerSequence,
				  bool bTopStrand,
				  Vector<int> &vSiteCounts,
				  Vector<int> &vSites)
{
  Vector<Sequence> vResult;
  int i;

  Motif hotMotif, coldMotif, n1Motif, n2Motif;
  if(bTopStrand) {

    Motif wrcMotif(3,2); wrcMotif[0] = 9; wrcMotif[1] = 5; wrcMotif[2] = C; // WRC
    Motif sycMotif(3,2); sycMotif[0] = 4; sycMotif[1] = 8; sycMotif[2] = C; // SYC
    Motif srcMotif(3,2); srcMotif[0] = 4; srcMotif[1] = 5; srcMotif[2] = C; // SRC (Neutral type 1)
    Motif wycMotif(3,2); wycMotif[0] = 9; wycMotif[1] = 8; wycMotif[2] = C; // WYC (Neutral type 2)

    hotMotif = wrcMotif;
    coldMotif = sycMotif;
    n1Motif = srcMotif;
    n2Motif = wycMotif;

  } else {

    Motif gywMotif(3,0); gywMotif[0] = G; gywMotif[1] = 8; gywMotif[2] = 9; // GYW
    Motif grsMotif(3,0); grsMotif[0] = G; grsMotif[1] = 5; grsMotif[2] = 4; // GRS
    Motif grwMotif(3,0); grwMotif[0] = G; grwMotif[1] = 5; grwMotif[2] = 9; // GRW (Neutral type 1)
    Motif gysMotif(3,0); gysMotif[0] = G; gysMotif[1] = 8; gysMotif[2] = 4; // GYS (Neutral type 2)

    hotMotif = gywMotif;
    coldMotif = grsMotif;
    n1Motif = gysMotif;
    n2Motif = grsMotif;

  }

  MotifIdentifier hotMotifIdentifier(hotMotif);
  vector<int> hotPositions = hotMotifIdentifier.shiftedFind(sd.vRef);

  MotifIdentifier coldMotifIdentifier(coldMotif);
  vector<int> coldPositions = coldMotifIdentifier.shiftedFind(sd.vRef);

  MotifIdentifier n1MotifIdentifier(n1Motif);
  vector<int> neutralPositions = n1MotifIdentifier.shiftedFind(sd.vRef);

  MotifIdentifier n2MotifIdentifier(n2Motif);
  vector<int> n2Positions = n2MotifIdentifier.shiftedFind(sd.vRef);

  neutralPositions.insert( neutralPositions.end(), n2Positions.begin(), n2Positions.end() );

  vector<int> allSites(neutralPositions);
  allSites.insert( allSites.end(), hotPositions.begin(), hotPositions.end() );
  allSites.insert( allSites.end(), coldPositions.begin(), coldPositions.end() );
  vSites.newsize(allSites.size());
  copy(allSites.begin(),allSites.end(),vSites.begin());

  double hotMutability = extractMeanMutability(vMutationProfile,hotPositions);
  double coldMutability = extractMeanMutability(vMutationProfile,coldPositions);
  double neutralMutability = extractMeanMutability(vMutationProfile,neutralPositions);

  Vector<double> vResults(3);
  vResults[HOT] = hotMutability;
  vResults[COLD] = coldMutability;
  vResults[NEUTRAL] = neutralMutability;

  vSiteCounts.newsize(3);
  vSiteCounts[HOT] = hotPositions.size();
  vSiteCounts[COLD] = coldPositions.size();
  vSiteCounts[NEUTRAL] = neutralPositions.size();

  return(vResults);
}
