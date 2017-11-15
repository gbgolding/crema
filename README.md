# lncRNA prediction by machine learning classifiers

Classifying transcripts as lncRNAs is difficult because there is no consensus
on what a lncRNA truly is.
Instead of using thresholds, or rules, for identifying lncRNAs, this tool uses
an ensemble stacking method of 8 different gradient boostling models to
predict lncRNAs.
Trained only on true, validated lncRNAs, this method has been tested on plant
transcriptomes with a high success rate.

See out publication at:


## Getting started

To use this tool, simply clone this repository on your machine by:
```
git clone https://github.com/caitsimop/lncRNApredTool.git
```

### Prerequisites

To use this tool you will need:

1. python3 
    - [biopython](http://biopython.org/)
    - [scikit-learn](http://scikit-learn.org)  
2. [CPAT v1.2.1](http://rna-cpat.sourceforge.net/)
    - python2
3. [DIAMOND](https://github.com/bbuchfink/diamond)

## Example

Before you can run tool, you'll need to remove all rRNAs and tRNAs from your
input data. 

Then, you will need to run cpat.py. An example:

```
cpat.py -g your_transcript_fasta_file.fa -o cpat_output.txt -x ./cpat_models/ath_hexamer -d ./cpat_models/ath_logit.RData
```

Then, DIAMOND:

```
diamond blastx -d swissprot.dmnd -q your_transcript_fasta_file.fa -o diamond_output.txt - e 0.001 -k 5 --matrix BLOSUM62 --gapopen 11 --gapextend 1 --more-sensitive -f 6 qseqid pident length qframe qstart qend sstart send evalue bitscore
```

Once you have identified your transcript features using CPAT and DIAMOND, you
can run the tool!

```
python3 bin/predict.py -f your_transcript_fasta_file.fa -c cpat_output.txt -d diamond_output.txt
```
