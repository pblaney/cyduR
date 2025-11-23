#include <cmath>

#include "scharff_utils.h"
#include "analysis.h"
#include "DiscreteHistogram.h"
#include "Sequence.h"
#include "MotifIdentifier.h"
#include "DiscreteHistogram.h"
#include "MotifMutationPair.h"
#include "StatsSampler.h"
#include "RandomizedIota.h"
#include "SimpleFastaReader.h"
#include "EgnetProperties.h"

void 
generateMutationHistograms(map<string,string> mss, 
			   string sRef, 
			   int maxPosition, 
			   int cdrCutoff,
			   string sOutputDir)
{
  int i,j;
  createDirectory(sOutputDir);
  Sequence vRef = convertSequence(sRef);
  Vector< Sequence > vvSeq = convertSequence(mss);
  chopSequences(vRef,vvSeq,maxPosition);
  //removeGapsInReferenceSequence(vRef,vvSeq);
  prechopSequences(vRef,vvSeq,cdrCutoff);
  maxPosition = vRef.size();

  char mutations[] = {'G','C','A','T'};
  int mutation;
  vector<int> vCounter(4,0);

  FOR(mutation,4) {
    char mutationFrom = mutations[mutation];
    string sFileName("from.X.dat");
    sFileName[5] = mutationFrom;
    ofstream outFile((sOutputDir+sFileName).c_str());
    FOR(i,maxPosition) {
      DiscreteHistogram dh;
      if( vRef[i]==mutation ) {
	FOR(j,vvSeq.size()) {
	  if( (vRef[i]!=vvSeq[j][i]) && (vRef[i]!=GAP) && (vvSeq[j][i]!=GAP) ) { // mutation at position i
	    dh.sample(vvSeq[j][i]);
	    vCounter[mutation]++;
	  }
	}
      }
      outFile << TB(i);
      FOR(j,4) {
	outFile << TB((double)dh[j]/vvSeq.size());
      }
      outFile << endl;
    }
    outFile.close();

  }

  FOR(mutation,4)    {
    char mutationFrom = mutations[mutation];
    string sFileName("from.X.total.dat");
    sFileName[5] = mutationFrom;
    ofstream outFile((sOutputDir+sFileName).c_str());
    outFile << vCounter[mutation] << endl;
    outFile.close();
  }

  Motif wrcMotif(3,0,2); wrcMotif[0] = 9; wrcMotif[1] = 5; wrcMotif[2] = 1; // WRC
  MotifIdentifier wrcMotifIdentifier(wrcMotif);
  vector<int> vMotifPositions   = wrcMotifIdentifier.shiftedFind(vRef);
  ofstream outFile((sOutputDir+string("wrc.motif.positions.dat")).c_str());
  FOR(i,vMotifPositions.size())
    outFile << vMotifPositions[i] << endl;
  outFile.close();

  Motif gywMotif(3,0,0); gywMotif[0] = 0; gywMotif[1] = 8; gywMotif[2] = 9; // GYW
  MotifIdentifier gywMotifIdentifier(gywMotif);
  vMotifPositions = gywMotifIdentifier.shiftedFind(vRef);
  outFile.open((sOutputDir+string("gyw.motif.positions.dat")).c_str());
  FOR(i,vMotifPositions.size())
    outFile << vMotifPositions[i] << endl;
  outFile.close();

  {
    Motif motif(1,0,0); motif[0] = 1; // C
    MotifIdentifier motifIdentifier(motif);
    vector<int> vMotifPositions   = motifIdentifier.shiftedFind(vRef);
    ofstream outFile((sOutputDir+string("c.motif.positions.dat")).c_str());
    FOR(i,vMotifPositions.size())
      outFile << vMotifPositions[i] << endl;
    outFile.close();
  }

  {
    Motif motif(1,0,0); motif[0] = 0; // G
    MotifIdentifier motifIdentifier(motif);
    vector<int> vMotifPositions   = motifIdentifier.shiftedFind(vRef);
    ofstream outFile((sOutputDir+string("g.motif.positions.dat")).c_str());
    FOR(i,vMotifPositions.size())
      outFile << vMotifPositions[i] << endl;
    outFile.close();
  }

  {
    Motif motif(2,0,1); motif[0]  = 9; motif[1]  = 2; 
    MotifIdentifier motifIdentifier(motif);
    vector<int> vMotifPositions   = motifIdentifier.shiftedFind(vRef);
    ofstream outFile((sOutputDir+string("wa.motif.positions.dat")).c_str());
    FOR(i,vMotifPositions.size())
      outFile << vMotifPositions[i] << endl;
    outFile.close();
  }

  {
    Motif motif(2,0,0); motif[1]  = 2; motif[1]  = 9;  
    MotifIdentifier motifIdentifier(motif);
    vector<int> vMotifPositions   = motifIdentifier.shiftedFind(vRef);
    ofstream outFile((sOutputDir+string("tw.motif.positions.dat")).c_str());
    FOR(i,vMotifPositions.size())
      outFile << vMotifPositions[i] << endl;
    outFile.close();
  }

  { // also write the number of mutations per mutated sequence
    ofstream outFile((sOutputDir+string("mutations.per.sequence.dat")).c_str());
    FOR(i,vvSeq.size())
      outFile << vvSeq[i].distance(vRef) << endl;
    outFile.close();
  }

}

Matrix<DiscreteHistogram>
generateInternalMutationHistograms(Sequence vRef,
				   Vector< Sequence > vvSeq)
{
  int i,j;
  int maxPosition = vRef.size();
  int mutation;
  Matrix<DiscreteHistogram> vResult(5,maxPosition);

  FOR(mutation,4) {
    FOR(i,maxPosition) {
      if( vRef[i]==mutation ) {
	FOR(j,vvSeq.size()) {
	  if( vRef[i] != vvSeq[j][i] ) // mutation at position i
	    vResult[mutation][i].sample(vvSeq[j][i]);
	}
      }
    }
  }
  return vResult;
}

vector< Matrix<int> >
generateInternalTransitionMatrices(Sequence vRef,
				   Vector< Sequence > vvSeq)
{
  int i,j;
  int maxPosition = vRef.size();
  int mutation;
  vector< Matrix<int> > vResult(maxPosition, Matrix<int>(5,5,0) );

  FOR(i,maxPosition) {
    FOR(j,vvSeq.size()) {
      int from = vRef[i];
      int to = vvSeq[j][i];
      if( from != to ) // mutation at position i
	vResult[i][from][to]++;
    }
  }
  return vResult;
}

vector< Matrix<int> >
generateInternalTransitionMatrices(SequenceDataset sequenceDataset)
{
  return generateInternalTransitionMatrices(sequenceDataset.vRef,sequenceDataset.vvSeq);
}

vector< Matrix<double> >
generateNormalizedTransitionMatrices(Sequence vRef,
				     Vector< Sequence > vvSeq)
{
  int i;
  vector< Matrix<int> > r1 = generateInternalTransitionMatrices(vRef,vvSeq);
  vector< Matrix<double> > r2(r1.size());
  FOR(i,r1.size())
    r2[i] = normalize(r1[i]);
  return r2;
}

double 
compareTransitionMatrices(const vector< Matrix<double> > &a,
			  const vector< Matrix<double> > &b)
{
  int i,j,k;
  double diff,r=0;
  FOR(i,a.size()) {
    FOR(j,a[i].num_rows()) {
      FOR(k,a[i].num_cols()) {
	diff = a[i][j][k]-b[i][j][k];
	r += diff*diff;
      }
    }
  }
  r = r/a.size();
  return r;
}

vector< Matrix<double> >
generateDetailedProfile(Sequence vRef,
			Vector< Sequence > vvSeq)
{
  int i,j,k,sum=0;
  vector< Matrix<int> > r1 = generateInternalTransitionMatrices(vRef,vvSeq);
  vector< Matrix<double> > r2(r1.size(),Matrix<double>(5,5,0.0));
  /*
  FOR(i,r1.size()) {
    FOR(j,r1[i].num_rows()) {
      FOR(k,r1[i].num_cols()) {
	sum += r1[i][j][k];
      }
    }
  }
  */
  FOR(i,r1.size()) {
    FOR(j,r1[i].num_rows()) {
      FOR(k,r1[i].num_cols()) {
	r2[i][j][k] = (double)r1[i][j][k] / vvSeq.size();
      }
    }
  }
  return r2;
}

//enum {GC,GA,GT,CG,CA,CT,AG,AC,AT,TG,TC,TA};

class MutationClassification
{
public:
  int position, from, to , type;
  MutationClassification(int _position, int _from, int _to)
    : position(_position), from(_from), to(_to)
  {
  }
  ~MutationClassification()
  {
  }
  int dist(const MutationClassification &m) {
    return abs(m.position-position);
  }
};

enum {LEFT,RIGHT};

class RlPositions : public vector< set<int> >
{
public:
  RlPositions()
    : vector< set<int> >(2)
  {}
  ~RlPositions()
  {}
};

class MutationMap : public map<int,RlPositions>
{
public:
  void insertAT(int x) {
    // 1. find key closest to x
    MutationMap::iterator it, closest = this->begin();
    map<int,RlPositions> &_this = *this;

    IFOR(it,_this) {
      int keyPos = it->first;
      //PR(keyPos);
      int d = keyPos-x;
      int currentClosest = (closest->first)-x;
      if( abs(d) < abs(currentClosest) )
	closest = it;
    }

    // 2. insert x
    RlPositions &rlPositions = closest->second;
    int d1 = x - closest->first;

    if(d1<0)
      rlPositions[LEFT].insert(-d1);
    else 
      rlPositions[RIGHT].insert(d1);

    if( ++closest!=_this.end() ) {
      int d2 = x - closest->first;
      //cerr << "d1=" << d1 << "\td2=" << d2 << endl;
      if((-d2)==d1) {
	RlPositions &rlPositions = closest->second;
	assert(d2<0);
	rlPositions[LEFT].insert(-d2);
      }
    }


  }
  
};

void 
generateNeighbouringHistograms(map<string,string> mss, string sRef, 
			       int maxPosition, 
			       string sOutputDir, int cdrCutoff)
{
  int i,j,k;

  Sequence vRef = convertSequence(sRef);
  Vector< Sequence > vvSeq = convertSequence(mss);
  chopSequences(vRef,vvSeq,maxPosition);
  removeGapsInReferenceSequence(vRef,vvSeq);
  prechopSequences(vRef,vvSeq,cdrCutoff);
  maxPosition = vRef.size();
  DiscreteHistogram lHist, rHist;

  FOR(i,vvSeq.size()) {
    
    MutationMap mutationMap;
    Sequence mutSeq = vvSeq[i];

    FOR(j,maxPosition) {
      //if( (vRef[j] != mutSeq[j]) && ( (vRef[j] == C) || (vRef[j] == G) ) )
      if( ((vRef[j] == C) && (mutSeq[j] == T)) || ((vRef[j] == G) && (mutSeq[j] == A)) )
	mutationMap[j] = RlPositions();
    }

    if(mutationMap.size()) {
      FOR(j,maxPosition) {
	if( (vRef[j] != mutSeq[j]) && ( (vRef[j] == A) || (vRef[j] == T) ) )
	  mutationMap.insertAT(j);
      }
    }

    MutationMap::const_iterator it;
    set<int>::const_iterator itI;
    IFOR(it,mutationMap) {
      const RlPositions &rlPositions = it->second;

      IFOR(itI,rlPositions[LEFT])
	lHist.sample( *itI );
      IFOR(itI,rlPositions[RIGHT])
	rHist.sample( *itI );
      
      /*
      if(rlPositions[LEFT].size()) 
	lHist.sample( *(rlPositions[LEFT].begin()) );
      if(rlPositions[RIGHT].size())
	rHist.sample( *(rlPositions[RIGHT].begin()) );
      */
    }

  }

  createDirectory(sOutputDir);
  ofstream outFile( (sOutputDir+"out.dat").c_str() );
  for(i=lHist.getUpperBound()-1; i>=1; i--)
    outFile << TB(-i) << lHist[i] << endl;
  for(i=1; i<rHist.getUpperBound(); i++)
    outFile << TB(i) << rHist[i] << endl;
  outFile.close();

}


void 
generateCDRGroupHistograms(map<string,string> mss, string sRef, 
			   int maxPosition, 
			   string sOutputDir, int cdrCutoff)
{
  int i,j,k;

  Sequence _vRef = convertSequence(sRef);
  Vector< Sequence > vvSeq = convertSequence(mss);
  removeGapsInReferenceSequence(_vRef,vvSeq);
  chopSequences(_vRef,vvSeq,maxPosition);
  DiscreteHistogram lHist, rHist;
  int rseqLength = maxPosition-cdrCutoff;
  map<Sequence, map<Sequence,int> > mCdr;

  PR(rseqLength);

  Sequence vRef(rseqLength,0);
  FOR(i,rseqLength)
    vRef[i] = _vRef[cdrCutoff+i];

  FOR(i,vvSeq.size()) {
    // 1. cut key
    Sequence k(cdrCutoff,0), v(rseqLength,0);
    FOR(j,maxPosition) {
      if(j < cdrCutoff)
	k[j] = vvSeq[i][j];
      else
	v[j-cdrCutoff] = vvSeq[i][j];
    }
    
    mCdr[k][v] += 1;
  }

  map<Sequence, map<Sequence,int> >::iterator itA;
  map<Sequence,int>::iterator itB;

  createDirectory(sOutputDir);
  ofstream outFile( (sOutputDir+"out.dat").c_str() );

  PR(mCdr.size());

  IFOR(itA,mCdr) {

    DiscreteHistogram lHist, rHist;

    IFOR(itB,itA->second) {

      PR(itB->second);

      MutationMap mutationMap;
      Sequence mutSeq = itB->first;

      FOR(j,rseqLength) {
	if( (vRef[j] != mutSeq[j]) && ( (vRef[j] == C) || (vRef[j] == G) ) )
	  mutationMap[j] = RlPositions();
      }

      if(mutationMap.size()) {
	FOR(j,rseqLength) {
	  if( (vRef[j] != mutSeq[j]) && ( (vRef[j] == A) || (vRef[j] == T) ) )
	    mutationMap.insertAT(j);
	}
      }

      MutationMap::const_iterator it;
      set<int>::const_iterator itI;
      IFOR(it,mutationMap) {
	const RlPositions &rlPositions = it->second;
	
	IFOR(itI,rlPositions[LEFT])
	  lHist.sample( *itI );
	IFOR(itI,rlPositions[RIGHT])
	  rHist.sample( *itI );
      }
      
    }

    for(i=rseqLength; i>=1; i--)
      outFile << TB(lHist[i]) ;
    FOR(i,rseqLength)
      outFile << TB(rHist[i]) ;
    outFile << endl;

    cout << endl;

    
  }

  outFile.close();
}


void 
generateHistogramsUsingMSH2Hotspots(map<string,string> mss, string sRef, 
				    int maxPosition, 
				    string sOutputDir, int cdrCutoff)
{
  int i,j,k;

  Sequence vRef = convertSequence(sRef);
  map<string,string> mssMsh, copyMss(mss);
  map<string,string>::iterator itI;
  
  IFOR(itI,copyMss) {
    string sK = itI->first;
    if(sK.find(">UNGMSH") == 0) {
      mssMsh.insert(*itI);
      mss.erase(sK);
    }
  }

  Vector< Sequence > vvSeqMsh = convertSequence(mssMsh);

  Sequence copy_vRef(vRef);
  prechopSequences(copy_vRef,vvSeqMsh,cdrCutoff);
  chopSequences(copy_vRef,vvSeqMsh,maxPosition-cdrCutoff);

  // we now want to extract the mutation positions
  set<int> vMshMutationPositions;

  FOR(i,vvSeqMsh.size()) {
    FOR(j,vvSeqMsh[i].size())
      if(copy_vRef[j] != vvSeqMsh[i][j]) {
	vMshMutationPositions.insert(j);
      }
  }

  Vector< Sequence > vvSeq = convertSequence(mss);
  prechopSequences(vRef,vvSeq,cdrCutoff);
  chopSequences(vRef,vvSeq,maxPosition-cdrCutoff);  
  maxPosition = vRef.size();

  DiscreteHistogram lHist, rHist;

  FOR(i,vvSeq.size()) {
    
    MutationMap mutationMap;
    Sequence mutSeq = vvSeq[i];
    set<int>::iterator itJ;

    FOR(j,maxPosition) {
      //if( ((vRef[j] == C) && (mutSeq[j] != vRef[j])) || ((vRef[j] == G) && (mutSeq[j] == A)) && 
      if( (vRef[j] != mutSeq[j]) && ( (vRef[j] == C) || (vRef[j] == G) )  && 
	  (vMshMutationPositions.find(j)!=vMshMutationPositions.end()) )
	mutationMap[j] = RlPositions();
    }

    if(mutationMap.size()) {
      FOR(j,maxPosition) {
	if( (vRef[j] != mutSeq[j]) && ( (vRef[j] == A) || (vRef[j] == T) ) )
	  mutationMap.insertAT(j);
      }
    }

    MutationMap::const_iterator it;
    set<int>::const_iterator itI;
    IFOR(it,mutationMap) {
      const RlPositions &rlPositions = it->second;

      IFOR(itI,rlPositions[LEFT])
	lHist.sample( *itI );
      IFOR(itI,rlPositions[RIGHT])
	rHist.sample( *itI );
      
    }

  }

  createDirectory(sOutputDir);
  ofstream outFile( (sOutputDir+"out.dat").c_str() );
  for(i=lHist.getUpperBound()-1; i>=1; i--)
    outFile << TB(-i) << lHist[i] << endl;
  for(i=1; i<rHist.getUpperBound(); i++)
    outFile << TB(i) << rHist[i] << endl;
  outFile.close();

}


void 
testProximityHypothesis(map<string,string> mss, 
			string sRef, 
			int maxPosition, 
			int cdrCutoff,
			string sOutputDir)
{
  int i,j;

  Sequence vRef = convertSequence(sRef);
  Vector< Sequence > vvSeq = convertSequence(mss);

  chopSequences(vRef,vvSeq,maxPosition);
  //removeGapsInReferenceSequence(vRef,vvSeq);
  prechopSequences(vRef,vvSeq,cdrCutoff);

  set<int> mutationPositions;
  set<int>::iterator it;
  maxPosition = vRef.size();
  char mutations[] = {'G','C','A','T'};

  /*
  FOR(i,vvSeq.size()) {    
    Sequence mutSeq = vvSeq[i];
    FOR(j,mutSeq.size()) {
      if(vRef[j]==G && mutSeq[j]==A)
	mutationPositions.insert(j);
    }
  }

  vector<DiscreteHistogram> vHist(5);
  
  IFOR(it,mutationPositions) {
    int pos = *it;
    if(pos>1 && pos <(maxPosition-2)) {
      int initPos = pos-2;
      FOR(i,5)
	vHist[i].sample(vRef[initPos+i]);
    }
  }

  FOR(i,5) {
    cerr << "position " << i << endl; 
    FOR(j,4) {
      cerr << TB((int)vHist[i][j]);
    }
    cerr << endl;
  }
  exit(0);
  */

  int mutation;

  //Motif wrcMotif(3,0,2); wrcMotif[0] = 9; wrcMotif[1] = 5; wrcMotif[2] = 1; // WRC
  Motif wrcMotif(3,0,2); wrcMotif[0] = 9; wrcMotif[1] = 0; wrcMotif[2] = 1; // WGC
  //Motif wrcMotif(2,0,1); wrcMotif[0] = 0; wrcMotif[1] = 1; // GC
  MotifIdentifier motifIdentifier(wrcMotif);

  vector<int> vMotifPositions   = motifIdentifier.find(vRef);

  //PR(vMotifPositions.size());

  int numSites = 0, cMutations=0, gMutations=0, bothMutations=0;
  int total_numSites = 0, total_cMutations=0, total_gMutations=0, total_bothMutations=0;

  FOR(j,vMotifPositions.size()) {

    int gPosition = vMotifPositions[j]+1;
    int cPosition = vMotifPositions[j]+2;

    FOR(i,vvSeq.size()) {

      Sequence mutSeq = vvSeq[i];
      if( (mutSeq[gPosition] != GAP) && (mutSeq[cPosition] != GAP) ) {

	numSites++;

	if(mutSeq[cPosition] != C)
	  cMutations++;

	if(mutSeq[gPosition] != G)
	  gMutations++;

	if( (mutSeq[cPosition] != C) && (mutSeq[gPosition] != G) ) // occur together in the same mutant sequence
	  bothMutations++;
      
      }
    }

    total_numSites += numSites;
    total_cMutations += cMutations; 
    total_gMutations += gMutations;
    total_bothMutations += bothMutations;

    double cFrequency = (double)cMutations/numSites;
    double gFrequency = (double)gMutations/numSites;
    double bothFrequency = (double)bothMutations/numSites;
    cerr << TB(cdrCutoff+vMotifPositions[j]) << TB(cFrequency) << TB(gFrequency) 
	 << TB(bothFrequency) << TB(gFrequency*cFrequency) << endl;
  }

  double total_cFrequency = (double)total_cMutations/total_numSites;
  double total_gFrequency = (double)total_gMutations/total_numSites;
  double total_bothFrequency = (double)total_bothMutations/total_numSites;
  cerr << "Global\t\n\t" << TB(total_cFrequency) << TB(total_gFrequency) 
    //<< TB(total_bothMutations)
       << TB(total_bothFrequency) << TB(total_gFrequency*total_cFrequency) << endl;
  
  int estimatedProduct = (int)round( (total_gFrequency*total_cFrequency) * total_numSites );
  cerr << "chisq\t" << total_bothMutations << "\t" << estimatedProduct << "\t" << total_numSites << endl;
}


void 
analyzeMotifAdvantage(map<string,string> mss, 
		      string sRef, 
		      int maxPosition, 
		      int cdrCutoff,
		      string sOutputDir)
{
  int i,j;

  Sequence vRef = convertSequence(sRef);
  Vector< Sequence > vvSeq = convertSequence(mss);

  chopSequences(vRef,vvSeq,maxPosition);
  prechopSequences(vRef,vvSeq,cdrCutoff);

  createDirectory(sOutputDir);

  Vector<double> vSmp = buildSiteMutationProfile(vRef, vvSeq);

  {
    Motif motif(1,0,0); motif[0] = 1; // C
    MotifIdentifier motifIdentifier(motif);
    vector<int> vSites = motifIdentifier.shiftedFind(vRef);
    ofstream outFile( (sOutputDir+"C.dat").c_str() );
    FOR(i,vSites.size())
      outFile << vSmp[ vSites[i] ] << endl;
    outFile.close();
  }

  {
    Motif motif(1,0,0); motif[0] = 0; // G
    MotifIdentifier motifIdentifier(motif);
    vector<int> vSites = motifIdentifier.shiftedFind(vRef);
    ofstream outFile( (sOutputDir+"G.dat").c_str() );
    FOR(i,vSites.size())
      outFile << vSmp[ vSites[i] ] << endl;
    outFile.close();
  }

  {
    Motif motif(3,0,2); motif[0] = 9; motif[1] = 5; motif[2] = 1; // WRC
    MotifIdentifier motifIdentifier(motif);
    vector<int> vSites = motifIdentifier.shiftedFind(vRef);
    ofstream outFile( (sOutputDir+"WRC.dat").c_str() );
    FOR(i,vSites.size())
      outFile << vSmp[ vSites[i] ] << endl;
    outFile.close();
  }

  {
    Motif motif(3,0,0); motif[0] = 0; motif[1] = 8; motif[2] = 9; // GYW
    MotifIdentifier motifIdentifier(motif);
    vector<int> vSites = motifIdentifier.shiftedFind(vRef);
    ofstream outFile( (sOutputDir+"GYW.dat").c_str() );
    FOR(i,vSites.size())
      outFile << vSmp[ vSites[i] ] << endl;
    outFile.close();
  }

}

vector<int>
countMotifs(const vector<MotifMutationPair> &vMotifMutationPair, 
	    const Sequence &vAux)
{
  vector<int> vResult(vMotifMutationPair.size());
  int pairChoice;
  FOR(pairChoice,vMotifMutationPair.size()) { // need to exclude the Null mutation
    Motif curMotif = vMotifMutationPair[pairChoice].motif;
    MotifIdentifier motifIdentifier(curMotif);
    vector<int> vMotifPositions = motifIdentifier.shiftedFind(vAux);
    //cerr << TB(pairChoice) << vMotifPositions.size() << endl;
    vResult[pairChoice] = vMotifPositions.size();
  }
  return vResult;
}

void
analyzeMotifDistInRandomizedSeqs(map<string,string> mss, string sRef, 
				 string sOutputDir, int nodeID, int maxPosition, int cdrCutoff,
				 int numRandomizations)
{
  int i,j,k,l,m;
  Sequence vRef = convertSequence(sRef);
  Vector< Sequence > vvSeq = convertSequence(mss);
  chopSequences(vRef,vvSeq,maxPosition);
  prechopSequences(vRef,vvSeq,cdrCutoff);

  vector<MotifMutationPair> vMotifMutationPair;

  // WARNING: Do not change the order of these, or BubbleEvolver won't work properly
  {
    Motif motif(1,0,0); motif[0] = 1; // C
    vMotifMutationPair.push_back( MotifMutationPair(motif,createMutation( Vector<char>(2,"CT") )) );
  }

  {
    Motif motif(3,0,2); motif[0] = 9; motif[1] = 5; motif[2] = 1; // WRC
    vMotifMutationPair.push_back( MotifMutationPair(motif,createMutation( Vector<char>(2,"CT") )) );
  }

  {
    Motif motif(1,0,0); motif[0] = 0; // G
    vMotifMutationPair.push_back( MotifMutationPair(motif,createMutation( Vector<char>(2,"GA") )) );
  }

  {
    Motif motif(3,0,0); motif[0] = 0; motif[1] = 8; motif[2] = 9; // GYW
    vMotifMutationPair.push_back( MotifMutationPair(motif,createMutation( Vector<char>(2,"GA") )) );
  }

  vector<int> vResult = countMotifs(vMotifMutationPair,vRef);
  //cerr << vResult << endl;
  int actualWrcCount = vResult[1];
  int actualGywCount = vResult[3];
  RandomizedIota ri(vRef.size());
  int numTrials = 1000, wrcGTCount=0, gywGTCount=0;

  FOR(i,numTrials) {
    Sequence vAux(vRef);
    FOR(j,vAux.size()) 
      vAux[j] = vRef[ ri[j] ];
    vResult = countMotifs(vMotifMutationPair,vAux);
    int simWrcCount = vResult[1];
    if(actualWrcCount>simWrcCount)
      wrcGTCount++;
    int simGywCount = vResult[3];
    if(actualGywCount>simGywCount)
      gywGTCount++;
    ri.shuffle();
  }

  double wrcGTProp = (double)wrcGTCount/numTrials;
  double gywGTProp = (double)gywGTCount/numTrials;

  PR(wrcGTProp);
  PR(gywGTProp);
}

void
preprocessSequences(EgnetProperties propsConfig)
{
  int i,j;
  map<string,string> mssData, mssRefs;
  string sRefSequence;
  vector<SequenceDataset> vSequenceDataset;

  int rCutoff = atoi(propsConfig.getProperty("rCutoff","-1").c_str());
  string sInputFile = propsConfig.getProperty("inputFile");
  string sOutputFile = propsConfig.getProperty("outputFile");
  SimpleFastaReader sfr;
  mssData = sfr.read(sInputFile.c_str());
 
  vector<string> vRefKeys;
  EgnetProperties::iterator itAux = propsConfig.find("refKeys");
  while( (itAux!=propsConfig.end()) && (itAux->first=="refKeys"))
    vRefKeys.push_back( (itAux++)->second );

  Vector< Sequence > vvRefs(vRefKeys.size());
  FOR(i,vRefKeys.size()) {
    vvRefs[i] = convertSequence( mssData[vRefKeys[i]] );
    mssData.erase(vRefKeys[i]);
  }

  Vector< Sequence > vvData = convertSequence(mssData);

  chopSequences(vvData,rCutoff);
  chopSequences(vvRefs,rCutoff);
  
  Sequence sequenceMainRef = vvRefs[0];

  if(vvRefs.size() == 2) {
    Sequence sequenceAltRef  = vvRefs[1];  

    assert(sequenceMainRef.size()==sequenceAltRef.size());

    // 1. remove polymorphisms, so they aren't interpreted as mutations
    FOR(i,sequenceMainRef.size()) {
      if(sequenceMainRef[i] != sequenceAltRef[i]) { // polymorphism
	int mainBase = sequenceMainRef[i];
	int altBase = sequenceAltRef[i];
	FOR(j,vvData.size()) {
	  if(vvData[j][i] == altBase) // switch back to main base
	    vvData[j][i] = mainBase;
	}
      }
    }

  }

  removeGapsInReferenceSequence(sequenceMainRef,vvData);

  // 2. make any non-base entries in main data set same as mainref
  FOR(i,vvData.size()) {
    FOR(j,vvData[i].size()) {
      if(vvData[i][j]>3) { // i.e. not GCAT
	vvData[i][j] = sequenceMainRef[j];
      }
    }
  }

  removeUnmutatedSequences(sequenceMainRef,vvData);

  writeFastaSequences(sequenceMainRef, vvData, sOutputFile);

}
