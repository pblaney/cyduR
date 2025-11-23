#include <algorithm>
#include <sstream>
#include <gsl/gsl_cdf.h>

#include "HotspotStatisticsReporter.h"
#include "scharff_utils.h"
#include "BivariateNormalConditional.h"
#include "SimpleFastaReader.h"
#include "LineReader.h"

string
CombinedMotif::
getPrettyString()
{
  string r = bSenseStrandOnly ? "S_" : "";
  r += fwdMI.motif.getPrettyString();
  return(r);
}

void
HotspotStatisticsReporter::
addSearchMotif(string sMotif, int mutatedPosition, int bSenseStrand)
{
  int i, j, motifLength = sMotif.size();
  CombinedMotif combinedMotif;
  combinedMotif.bSenseStrandOnly = bSenseStrand;
  combinedMotif.label = sMotif;
  if(bSenseStrand)
    combinedMotif.label = "S_"+sMotif;

  Motif fwdMotif(motifLength,mutatedPosition);
  FOR(i,motifLength) {
    char thisLetter = sMotif[i];
    fwdMotif[i] = motifReference.mBP[thisLetter];
  }
  combinedMotif.fwdMI = MotifIdentifier(fwdMotif);

  if(!bSenseStrand) {  // It's both strands, so we need to do the other one too
    Motif revMotif(motifLength,motifLength-mutatedPosition-1);
    for(i=motifLength-1,j=0; i>=0; i--,j++) {
      char thisLetter = sMotif[i];
      int thisCode = motifReference.mBP[thisLetter];
      int revCode = motifReference.reverseCodes[thisCode];
      revMotif[j] = revCode;
    }
    combinedMotif.revMI = MotifIdentifier(revMotif);
  }

  vCM.push_back(combinedMotif);
}

HotspotStatisticsReporter::
HotspotStatisticsReporter()
  : bAssignedNucleotideProportions(false)
{
}

void
HotspotStatisticsReporter::
assignDefaultMotifs()
{
  addSearchMotif("WRC",2,BOTH_STRANDS);
  addSearchMotif("SYC",2,BOTH_STRANDS);
  //addSearchMotif("CCC",2,BOTH_STRANDS);
  //addSearchMotif("TC",1,BOTH_STRANDS);
  addSearchMotif("CG",0,SENSE_STRAND_ONLY);
  //addSearchMotif("TTC",2,BOTH_STRANDS);
  //addSearchMotif("TCC",2,BOTH_STRANDS);
  addSearchMotif("YC",1,BOTH_STRANDS);
  addSearchMotif("BC",1,BOTH_STRANDS);
  addSearchMotif("TA",0,SENSE_STRAND_ONLY);
  addSearchMotif("TT",0,SENSE_STRAND_ONLY);

  // RNA motifs
  addSearchMotif("ACA",1,SENSE_STRAND_ONLY); // APOBEC1 (Papavasiliou paper)
  addSearchMotif("CATC",3,SENSE_STRAND_ONLY); // APOBEC3/Monocytes paper (primary motif)
  addSearchMotif("CACC",3,SENSE_STRAND_ONLY); // APOBEC3/Monocytes paper (secondary)
  addSearchMotif("CCTC",3,SENSE_STRAND_ONLY); // APOBEC3/Monocytes paper (secondary)
  addSearchMotif("CTTC",3,SENSE_STRAND_ONLY); // APOBEC3/Monocytes paper (secondary)
  addSearchMotif("TATC",3,SENSE_STRAND_ONLY); // APOBEC3/Monocytes paper (secondary)

  //addSearchMotif("UAA",3,SENSE_STRAND_ONLY); // ADAR

  int i,j;
  FOR(i,4) { // add all NNC motifs
    FOR(j,4) {
      string sCandMotif;
      sCandMotif += MotifReference::bases[i];
      sCandMotif += MotifReference::bases[j];
      sCandMotif += "C";
      addSearchMotif(sCandMotif,2,BOTH_STRANDS);
    }
  }
  FOR(i,4) { // add all NC motifs
    string sCandMotif;
    sCandMotif += MotifReference::bases[i];
    sCandMotif += "C";
    addSearchMotif(sCandMotif,1,BOTH_STRANDS);
  }

  /*
  FOR(i,vCM.size()) // write out default motifs
    cout << vCM[i].fwdMI.motif.getPrettyString() << "\t" << (vCM[i].bSenseStrandOnly?"SENSE":"BOTH") << endl;
  exit(0);
  */
  
}

void
HotspotStatisticsReporter::
readExternalMotifs(ifstream &_motifStream)
{
  LineReader<string> motifLineReader(_motifStream);
  string sToken;
  char sDelim = '_';
  vector<string> vBuffer;
  while( motifLineReader.next(vBuffer) ) {
    string sMotif = vBuffer[0];
    int strand = (vBuffer[1] == "SENSE") ? SENSE_STRAND_ONLY : BOTH_STRANDS;
    
    istringstream tokenStream(sMotif);
    std::getline(tokenStream, sToken, sDelim);
    string sJoined(sToken);
    int motifPosition = sToken.length(); // position is length of 1st token
    while (std::getline(tokenStream, sToken, sDelim))
      sJoined += sToken;
    //cerr << TB(sJoined) << TB(motifPosition) << ( (strand==SENSE_STRAND_ONLY) ? "SENSE" : "BOTH" ) << endl;
    addSearchMotif(sJoined,motifPosition,strand);
  }
}

HotspotStatisticsReporter::
~HotspotStatisticsReporter()
{
}

char
HotspotStatisticsReporter::
getAA(const Sequence &sAux, int ntPosition)
{
  int codonStartPosition = ntPosition - (ntPosition%3);
  char cAa = '-';
  if( (codonStartPosition+3) <= sAux.size() ) {
    string sCodon;
    sCodon += motifReference.bases[ sAux[codonStartPosition] ];
    sCodon += motifReference.bases[ sAux[codonStartPosition+1] ];
    sCodon += motifReference.bases[ sAux[codonStartPosition+2] ];
    cAa = geneticCode.codonMap[sCodon];
  }
  return(cAa);
}

pair<char,string>
HotspotStatisticsReporter::
getAAwithCodon(const Sequence &sAux, int ntPosition)
{
  int codonStartPosition = ntPosition - (ntPosition%3);
  char cAa = '-';
  string sCodon;
  if( (codonStartPosition+3) <= sAux.size() ) {
    sCodon += motifReference.bases[ sAux[codonStartPosition] ];
    sCodon += motifReference.bases[ sAux[codonStartPosition+1] ];
    sCodon += motifReference.bases[ sAux[codonStartPosition+2] ];
    cAa = geneticCode.codonMap[sCodon];
  }
  return(make_pair(cAa,sCodon));
}

Sequence
HotspotStatisticsReporter::
randomizeCodons(const Sequence &sAux) 
{
  int i;
  // translate
  int numCodons = sAux.size()/3;
  string sAaSeq;
  FOR(i,numCodons)
    sAaSeq += getAA(sAux,i*3);

  string sNewSeq;
  FOR(i,sAaSeq.size()) {
    char thisAa = sAaSeq[i];
    vector<string> vCodons = geneticCode.aaMap[thisAa];
    //int randCodonChoice = getRandNum(0,vCodons.size()-1);
    FrequencyDependentRandomizer *auxFdr = mapFdr[thisAa];
    int randCodonChoice = auxFdr->getRand();
    string thisCodon = vCodons[randCodonChoice];
    sNewSeq += thisCodon;
  }

  Sequence sResult = convertSequence(sNewSeq);
  return(sResult);
}

Vector<int>
HotspotStatisticsReporter::
countHotspots(Sequence &sAux) 
{
  int i, numMotifs = vCM.size();
  Vector<int> vResult(numMotifs,0);
  FOR(i,numMotifs) {
    if(vCM[i].bSenseStrandOnly)
      vResult[i] = vCM[i].fwdMI.shiftedFind(sAux).size();  
    else {
      vResult[i] =  vCM[i].fwdMI.shiftedFind(sAux).size();  
      vResult[i] += vCM[i].revMI.shiftedFind(sAux).size();  
    }
  }
  return(vResult);
}


void 
HotspotStatisticsReporter::
countReplacements(Sequence &sAux, MotifIdentifier &mi, int mutationTo, int &numReplacements, int &numHotspots)
{
  int i;
  vector<int> vHotspots = mi.shiftedFind(sAux);
  numHotspots += vHotspots.size();

  FOR(i,vHotspots.size()) {
    int mutationPos = vHotspots[i];
    int mutationPosWithinCodon = mutationPos%3;
    pair<char,string> auxAC = getAAwithCodon(sAux,mutationPos);
    char prevAA = auxAC.first;
    string auxCodon = auxAC.second;

    auxCodon[mutationPosWithinCodon] = MotifReference::bases[mutationTo]; // auxCodon is now mutated
    char newAA = geneticCode.codonMap[auxCodon];

    if(newAA != prevAA)
      numReplacements++;
  }
  
}

Matrix<int>
HotspotStatisticsReporter::
countTransitionReplacements(Sequence &sAux) 
{
  int i, numMotifs = vCM.size();
  Matrix<int> mResult(numMotifs,2,0); // 1st column: # replacements, 2nd column: # hotspots

  FOR(i,numMotifs) {
    
    int frNuc = vCM[i].fwdMI.motif[vCM[i].fwdMI.motif.mutationPosition], toNuc, revNuc;
    switch(frNuc) {
    case C: toNuc = T; revNuc = A; break; // define all transitions, other strand also
    case T: toNuc = C; revNuc = G; break;
    case G: toNuc = A; revNuc = T; break;
    case A: toNuc = G; revNuc = C; break;
    default: fatalError("Unknown or degenerate \"from\" nucleotide");
    }
    
    if(vCM[i].bSenseStrandOnly)
      countReplacements(sAux,vCM[i].fwdMI,toNuc,mResult[i][0],mResult[i][1]);      
    else {
      countReplacements(sAux,vCM[i].fwdMI,toNuc,mResult[i][0],mResult[i][1]);
      countReplacements(sAux,vCM[i].revMI,revNuc,mResult[i][0],mResult[i][1]);
    }
  }

  return(mResult);
}


double
HotspotStatisticsReporter::
calcGcProportion(Sequence sOriginal) {
  int i, gcCount=0;
  FOR(i,sOriginal.size()) {
    if(sOriginal[i]==G || sOriginal[i]==C)
      gcCount++;
  }
  double r = (double)gcCount / sOriginal.size();
  return(r);
}

void
HotspotStatisticsReporter::
computeCodonWeights(double gcProportion) {

  int i,j;
  map<char,vector<string> >::iterator it;

  IFOR(it,geneticCode.aaMap) {

    char thisAA = it->first;
    vector<string> auxCodons = it->second;
    Vector<double> auxVf(auxCodons.size());

    FOR(i,auxCodons.size()) {
      string thisCodon = auxCodons[i];
      double v = 1.0;
      FOR(j,thisCodon.size()) {
	if( thisCodon[j]=='G' || thisCodon[j]=='C' )
	  v = v*gcProportion;
	else
	  v = v*(1.0-gcProportion);
      }
      auxVf[i] = v;
    }

    FrequencyDependentRandomizer *fdrPtr = new FrequencyDependentRandomizer(auxVf);
    fdrPtr->setSeed(0);
    mapFdr[thisAA] = fdrPtr;

  }

}

string 
HotspotStatisticsReporter::
changeNANtoNA(double x)
{
  string sResult;
  ostringstream aux;
  if(isnan(x))
    sResult = "NA";
  else {
    aux << x;
    sResult = aux.str();
  }
  return(sResult);
}

void 
HotspotStatisticsReporter:: 
evaluate(map<string,string> &mss, bool bCoding)
{
  int i,j,k;
  int numMotifs=vCM.size();
  map<string,string>::iterator it;
  Sequence sOriginal;
  vector<Sequence> vShuffled(mss.size()-1);

  i=0;
  IFOR(it,mss) {
    pair<string,string> pss = *it;
    string sKey = pss.first; // replicate#
    string sValue = pss.second; // shuffled sequence
    if ( sKey.find(">replicate") == 0 )
      vShuffled[i++] = convertSequence(sValue);
    else
      sOriginal = convertSequence(sValue);
  }
  
  Vector<int> originalNumHs = countHotspots(sOriginal);
  
  //cout << "original #:" << originalNumHs[2] << endl;
  
  Matrix<int> originalRp;
  Vector<double> originalRpFraction(numMotifs,0.0);

  if(bCoding) {
    originalRp = countTransitionReplacements(sOriginal);
    FOR(j,numMotifs) {
      if(originalRp[j][1] > 0)
	originalRpFraction[j] = (double) originalRp[j][0] / originalRp[j][1];
    }
  }

  //Vector<string> vMotifNames(numMotifs,"WRC SYC CCC TC CG TTC TCC YC BC S_TA S_TT");
  Vector<string> vMotifNames(numMotifs);
  FOR(j,numMotifs)
    vMotifNames[j] = vCM[j].label;

  Vector<int> vBelow(numMotifs,0), vBelowReplacement(numMotifs,0), vBelowReplacementFraction(numMotifs,0);
  Matrix<BivariateNormalConditional> mbnc(numMotifs,numMotifs), mRepTr(numMotifs,numMotifs), mRepTrFrac(numMotifs,numMotifs);
  //vss = vector of statistic samplers where each entry holds data for a different motif
  Vector<StatsSampler> vss(numMotifs), vssReplacement(numMotifs), vssReplacementFraction(numMotifs);
  Vector<DiscreteHistogram> vdh(numMotifs);
  Vector<int> adjustedNumIterations(numMotifs,vShuffled.size());

  FOR(i,vShuffled.size()) {
    Matrix<int> auxRp;

    Sequence sAux = vShuffled[i];
    if(bCoding)
      auxRp = countTransitionReplacements(sAux);

    Vector<int> auxNumHs = countHotspots(sAux);

    //cout << "simulated #:" << auxNumHs[2] << endl;

    
    FOR(j,numMotifs) {

      FOR(k,numMotifs) {
	mbnc[j][k].sample(auxNumHs[j],auxNumHs[k]);
      }

      vss[j].sample(auxNumHs[j]);
      if(auxNumHs[j] < originalNumHs[j])
	vBelow[j] += 1;

      if(bCoding) {

	vssReplacement[j].sample(auxRp[j][0]);
	if(auxRp[j][0] < originalRp[j][0])
	  vBelowReplacement[j] += 1;
	FOR(k,numMotifs) 
	  mRepTr[j][k].sample(auxRp[j][0],auxNumHs[k]);

	double auxFraction;   // NOTE THAT THE "ITERATIONS" USED FOR THE P-VALUE NEED TO BE ADJUSTED IF THE DENOMINATOR HERE IS ZERO
	if(auxRp[j][1] > 0) {
	  auxFraction = (double)auxRp[j][0] / auxRp[j][1];
	  vssReplacementFraction[j].sample(auxFraction);
	  if(auxFraction < originalRpFraction[j])
	    vBelowReplacementFraction[j] += 1;
	  FOR(k,numMotifs) 
	    mRepTrFrac[j][k].sample(auxFraction,auxNumHs[k]);
	} else
	  adjustedNumIterations[j]--;


      }
      //vdh[j].sample(auxNumHs[j]);

      //cerr << TB(auxNumHs[j]);
    }

    //cerr << endl;
  }
  //exit(0);

  bool bScreenOutput = false;
  double pBelow, pBelowApprox, observed, expected, expected_sd;

  FOR(j,numMotifs) {
    observed = originalNumHs[j];
    expected = vss[j].mean();
    expected_sd = vss[j].stddev();
    pBelow = (double) vBelow[j] / vShuffled.size();
    pBelowApprox = gsl_cdf_gaussian_P(observed-expected,expected_sd); //this is the p-value for the normally distributed data

    cout << "below"       << TB( vCM[j].getPrettyString() ) << changeNANtoNA(pBelow) << endl;
    cout << "belowApprox" << TB( vCM[j].getPrettyString() ) << changeNANtoNA(pBelowApprox) << endl;
    cout << "observed" << TB( vCM[j].getPrettyString() ) << changeNANtoNA(observed) << endl;
    cout << "expected" << TB( vCM[j].getPrettyString() ) << changeNANtoNA(expected) << endl;
    cout << "expectedSd" << TB( vCM[j].getPrettyString() ) << changeNANtoNA(expected_sd) << endl;
  }

  if(bCoding) {
    FOR(j,numMotifs) {
      observed = originalRp[j][0];
      expected = vssReplacement[j].mean();
      expected_sd = vssReplacement[j].stddev();
      pBelow = (double) vBelowReplacement[j] / vShuffled.size();
      pBelowApprox = gsl_cdf_gaussian_P(observed-expected,expected_sd);
      
      cout << "repTr_below" << TB( vCM[j].getPrettyString() ) << changeNANtoNA(pBelow) << endl;
      cout << "repTr_belowApprox" << TB( vCM[j].getPrettyString() ) << changeNANtoNA(pBelowApprox) << endl;
      cout << "repTr_observed" << TB( vCM[j].getPrettyString() ) << changeNANtoNA(observed) << endl;
      cout << "repTr_expected" << TB( vCM[j].getPrettyString() ) << changeNANtoNA(expected) << endl;
      cout << "repTr_expectedSd" << TB( vCM[j].getPrettyString() ) << changeNANtoNA(expected_sd) << endl;

      observed = originalRpFraction[j];
      expected = vssReplacementFraction[j].mean();
      expected_sd = vssReplacementFraction[j].stddev();
      pBelow = (double) vBelowReplacementFraction[j] / adjustedNumIterations[j];
      pBelowApprox = gsl_cdf_gaussian_P(observed-expected,expected_sd);
      
      cout << "repTrFrac_below" << TB( vCM[j].getPrettyString() ) << changeNANtoNA(pBelow) << endl;
      cout << "repTrFrac_belowApprox" << TB( vCM[j].getPrettyString() ) << changeNANtoNA(pBelowApprox) << endl;
      cout << "repTrFrac_observed" << TB( vCM[j].getPrettyString() ) << changeNANtoNA(observed) << endl;
      cout << "repTrFrac_expected" << TB( vCM[j].getPrettyString() ) << changeNANtoNA(expected) << endl;
      cout << "repTrFrac_expectedSd" << TB( vCM[j].getPrettyString() ) << changeNANtoNA(expected_sd) << endl;
    }
  }

  FOR(j,numMotifs) {
    FOR(k,j) {
      cout << "cor" << vCM[j].getPrettyString() << "x" << TB(vCM[k].getPrettyString()) << changeNANtoNA(mbnc[j][k].corrcoeff()) << endl;
      cout << "corRepTr" << vCM[j].getPrettyString() << "x" << TB(vCM[k].getPrettyString()) << changeNANtoNA(mRepTr[j][k].corrcoeff()) << endl;
      cout << "corRepTrFrac" << vCM[j].getPrettyString() << "x" << TB(vCM[k].getPrettyString()) << changeNANtoNA(mRepTrFrac[j][k].corrcoeff()) << endl;
    }
  }

  double pValue;
  FOR(j,numMotifs) {
    FOR(k,numMotifs) {
      if(j!=k) {
	pValue = mbnc[j][k].getV1ConditionalOnV2(originalNumHs[j],originalNumHs[k]);
	cout << "p" << vCM[j].getPrettyString() << "cond" << TB(vCM[k].getPrettyString()) << changeNANtoNA(pValue) << endl;
      }
      pValue = mRepTr[j][k].getV1ConditionalOnV2(originalRp[j][0],originalNumHs[k]);
      cout << "p" << vCM[j].getPrettyString() << "condRepTr" << TB(vCM[k].getPrettyString()) << changeNANtoNA(pValue) << endl;

      pValue = mRepTrFrac[j][k].getV1ConditionalOnV2(originalRpFraction[j],originalNumHs[k]);
      cout << "p" << vCM[j].getPrettyString() << "condRepTrFrac" << TB(vCM[k].getPrettyString()) << changeNANtoNA(pValue) << endl;

    }
  }

}

void 
HotspotStatisticsReporter::
assignNucleotideComposition(double aProportion, double cProportion, double gProportion)
{
  bAssignedNucleotideProportions = true;
  vProportionsACG = Vector<double>(3);
  vProportionsACG[0] = aProportion;
  vProportionsACG[1] = cProportion;
  vProportionsACG[2] = gProportion;
}

int main(int argc, char *argv[]) {
  HotspotStatisticsReporter hotspotStatisticsReporter;

  string sFilename(argv[1]);
  ifstream motifFile;
  if(argc > 2) {
    motifFile.open( argv[2] );
    if(!motifFile)
      fatalError("cannot open configuration file ...");
    hotspotStatisticsReporter.readExternalMotifs(motifFile);
  }
  else
    hotspotStatisticsReporter.assignDefaultMotifs();
    
  SimpleFastaReader sfr;
  map<string,string> mss = sfr.read(sFilename.c_str());
  hotspotStatisticsReporter.evaluate(mss,true);
}
