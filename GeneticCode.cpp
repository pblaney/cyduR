#include "GeneticCode.h"

GeneticCode::
GeneticCode()
{
  addEntry("GCT","A");
  addEntry("GCC","A");
  addEntry("GCA","A");
  addEntry("GCG","A");
  addEntry("TGT","C");
  addEntry("TGC","C");
  addEntry("GAT","D");
  addEntry("GAC","D");
  addEntry("GAA","E");
  addEntry("GAG","E");
  addEntry("TTT","F");
  addEntry("TTC","F");
  addEntry("GGT","G");
  addEntry("GGC","G");
  addEntry("GGA","G");
  addEntry("GGG","G");
  addEntry("CAT","H");
  addEntry("CAC","H");
  addEntry("ATT","I");
  addEntry("ATC","I");
  addEntry("ATA","I");
  addEntry("AAA","K");
  addEntry("AAG","K");
  addEntry("TTG","L");
  addEntry("TTA","L");
  addEntry("CTT","L");
  addEntry("CTC","L");
  addEntry("CTA","L");
  addEntry("CTG","L");
  addEntry("ATG","M");
  addEntry("AAT","N");
  addEntry("AAC","N");
  addEntry("CCT","P");
  addEntry("CCC","P");
  addEntry("CCA","P");
  addEntry("CCG","P");
  addEntry("CAA","Q");
  addEntry("CAG","Q");
  addEntry("CGT","R");
  addEntry("CGC","R");
  addEntry("CGA","R");
  addEntry("CGG","R");
  addEntry("AGA","R");
  addEntry("AGG","R");
  addEntry("TCT","S");
  addEntry("TCC","S");
  addEntry("TCA","S");
  addEntry("TCG","S");
  addEntry("AGT","S");
  addEntry("AGC","S");
  addEntry("ACT","T");
  addEntry("ACC","T");
  addEntry("ACA","T");
  addEntry("ACG","T");
  addEntry("GTT","V");
  addEntry("GTC","V");
  addEntry("GTA","V");
  addEntry("GTG","V");
  addEntry("TGG","W");
  addEntry("TAT","Y");
  addEntry("TAC","Y");
  addEntry("TAA","*");
  addEntry("TAG","*");
  addEntry("TGA","*");
}

GeneticCode::
~GeneticCode()
{
}

void
GeneticCode::
addEntry(char *_codon, char *_paa)
{
  string sCodon(_codon);
  char _aa = _paa[0];
  codonMap[sCodon] = _aa;
  map< char,vector<string> >::iterator it = aaMap.find(_aa);
  if(it == aaMap.end()) {
    vector<string> vAux(1,sCodon);
    aaMap[_aa] = vAux;
  } else {
    vector<string> &vAux = it->second;
    vAux.push_back(sCodon);
    //aaMap[_aa] = vAux;
  }
  
}
