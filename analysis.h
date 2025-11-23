#include "util.h"
#include "DiscreteHistogram.h"
#include "Sequence.h"
#include "SequenceDataset.h"

void 
generateMutationHistograms(map<string,string> mss, 
			   string sRef, 
			   int maxPosition, 
			   int cdrCutoff,
			   string sOutputDir);

Matrix<DiscreteHistogram>
generateInternalMutationHistograms(Sequence vRef,
				   Vector< Sequence > vvSeq);

void 
generateNeighbouringHistograms(map<string,string> mss, string sRef, 
			       int maxPosition, 
			       string sOutputDir, int cdrCutoff);

void 
generateCDRGroupHistograms(map<string,string> mss, string sRef, 
			   int maxPosition, 
			   string sOutputDir, int cdrCutoff);

void 
generateHistogramsUsingMSH2Hotspots(map<string,string> mss, string sRef, 
				    int maxPosition, 
				    string sOutputDir, int cdrCutoff);

void 
testProximityHypothesis(map<string,string> mss, 
			string sRef, 
			int maxPosition, 
			int cdrCutoff,
			string sOutputDir);

void 
analyzeMotifAdvantage(map<string,string> mss, 
		      string sRef, 
		      int maxPosition, 
		      int cdrCutoff,
		      string sOutputDir);

void
analyzeMotifDistInRandomizedSeqs(map<string,string> mss, string sRef, 
				 string sOutputDir, int nodeID, int maxPosition, int cdrCutoff,
				 int numRandomizations);

void
preprocessSequences(EgnetProperties propsConfig);

vector< Matrix<int> >
generateInternalTransitionMatrices(Sequence vRef,
				   Vector< Sequence > vvSeq);

vector< Matrix<int> >
generateInternalTransitionMatrices(SequenceDataset sequenceDataset);

vector< Matrix<double> >
generateNormalizedTransitionMatrices(Sequence vRef,
				     Vector< Sequence > vvSeq);

double 
compareTransitionMatrices(const vector< Matrix<double> > &a,
			  const vector< Matrix<double> > &b);

vector< Matrix<double> >
generateDetailedProfile(Sequence vRef,
			Vector< Sequence > vvSeq);
