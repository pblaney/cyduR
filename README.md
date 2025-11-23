# cytidine deaminase under-representation Reporter in R
>An updated tool-kit for over-/under- representation analysis of cytidine deaminase (AID and APOBEC) mutation motifs in a DNA sequence

### Introduction
The goal of `cyduR` is to extend the work from Maxwell Shapiro, Stephen Meier & Thomas MacCarthy that randomly shuffles a given FASTA sequence and then statistically tests for enrichment or depletion of cytidine deaminase NNC, WRC, SYC, and TC mutation motifs.

In `cyduR`, the shuffling algorithm and workflow orchestration scripts have been ported from `python` to `R` while maintaining the statistical reporter. The update primarily includes extension of the shuffling algorithm to any non-coding DNA sequence as well. Additional features for visualization have been included as well.

Citation: Shapiro, M., Meier, S., & MacCarthy, T. (2018). The cytidine deaminase under-representation reporter (CDUR) as a tool to study evolution of sequences under deaminase mutational pressure. *BMC Bioinformatics*, 19(1), 163. doi:10.1186/s12859-018-2161-y

Original Code Repository: [CDUR gitlab](https://gitlab.com/maccarthyslab/CDUR).




### Installation
CDUR has been updated to work on python3. It should also still work on python2. Using git, from your home directory, type:
```
git clone https://gitlab.com/maccarthyslab/CDUR.git
```
For ease of use, you should have the CDUR folder in your home directory. Otherwise, you may need to change where the reporter executable ("shmsim") is called.

1. Installation requires _pip_. To install _pip_ type:

```
sudo apt-get install python-pip
```

2. If you did not use `git clone`, download and unzip the CDUR repository and extract it (preferably create a CDUR folder - see above)
3. Run the _installCDUR.py_ file

```python
python installCDUR.py
```

### Usage
In order to run CDUR, one needs a fasta file with coding sequences (starting with ATG). 
If a fasta file contains multiple sequences, CDUR will analyze each one separately.

To run CDUR, type:

```python
python CDUR.py --arguments
```


The `--arguments` are as follows:

* `-i Input`
	- Specify the path to the fasta file you wish to analyze

* `-s Shuffling method choice`
    - The choices are gc3, n3, or dn23:
        * gc3 computes the GC content of the 3rd position in of each codon, then synonymously chooses a codon to maintain amino-acid integrity. This shuffle does not keep GC content, dinucleotide content, or codon bias well preserved
        
		* n3 shuffles all the 3rd positions codons and assigns then, whilst maintaining the underlying amino-acid sequence, assigns back each nulceotide randomly. This method maintains GC content, but may vary the dinucleotide frequency and codon pair bias.
        
		* dn23 computes the dinucleotide frequency of the sequence, and based on this frequency will synonymously choose 2nd and third positions for each codon. This does not maintain GC content, but will mostly preserve dinucleotide frequency and codon pair bias.
        
* `-r Number of shuffles`
    - Specifies the number of times the sequence is shuffled. Default r=1000.
	
* `-m Motif filename`
	- Specifies the file containing the mutation motifs to be analyzed together with  the strand (sense only or both strands). Within each motif the targeted nucleotide to analyze is delimited by the "\_", e.g., "AG\_C\_". Note that degenerate motifs are allowed, following the IUPAC nomenclature (e.g. W=A or T). The keywords "SENSE" or "BOTH" are used for strand. For example, to consider WRC motifs on both strands in the configuration file, one would require the line "WR\_C\_ &nbsp;&nbsp; BOTH" with a tab in between. Each motif to be considered should be on a new line. The motif configuration file "motif.txt" containing a default list of motifs is included.
    
* `-o Folder to place all output files`
	- Specify where to place the file containing the processed sequence, shuffled sequences, and results. If possible, use the full path.

* `-d Delete processed sequence and shuffled sequence`
	- Leave out if you want to keep all files. Otherwise, use `-d 1` to delete all but the results file.
	The "processed" sequence in this case may be from a multi-fasta file. Each file, for convenience, recieves its own fasta file, and then that file is opened and shuffled. This flag allows you to delete all of the intermediate files. 

* `-h Help`



### Example
Please go to <http://www.ams.sunysb.edu/~maccarth/software.html> and download AJ011405.1_protease.fas. This is the protease gene sequence for HIV-1.

```python
python CDUR.py -i AJ011405.1_protease.fas -s gc3 -r 2000
```

The output will be two files, both found in the directory from where CDUR is called:
* File 1. AJ011405.1_protease_gc3.fasta
	- This file contains the original sequence with all the shuffled sequences

* File 2. AJ011405.1_protease_gc3results.txt
	- This file contains a txt file with two columns. The first column contains the metric e.g. belowWRC, and the second column contains the P-value for under-representation, e.g. 0.639. The metrics are (all of the following consider forward and reverse compliment motifs, for example, counts of WRC and GYW):
		* below: % of shuffled sequences with fewer hotspots than original
		* observed: observed # hotspots in original sequence
		* expected: mean number of hotspots in shuffled sequences
		* expectedSd: standard deviation of hotspots in shuffled sequences
		* repTr: Number of non-synonymouse mutations
		* repTrFrac: (#non-synonymous mutations)/(#hotspots)
		* corXxY: correlation of X with Y
		* pXcondY: conditional of X given Y
                
In the calculation above, belowWRC yields the percentage of shuffled sequences that have fewer hotspots (WRC and GYW) than the original. As another example, in the example run above, repTR_belowWRC 0.503 means that the percentage of shuffled sequences with fewer non-synonymous mutations than the original is 0.503 (this constitutes an emprical P-value). Another example, corRepTrFracSYCxWRC 0.428467 means that the correlation of hotspots/non-syn between hotspots SYC and WRC is .0428467.
Note that if you run CDUR multiple times with the same filename, all the results are appended to the <sequence>results.txt file!

### Licenses
The CDUR program is licesned to Shapiro under the MIT license. The shuffling script was
licensed to danielmjorge copyright (c) 2015 under the MIT license and can be found at
<https://github.com/lauringlab/CodonShuffle>
