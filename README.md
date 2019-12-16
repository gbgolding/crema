# Classifying RNAs by Ensemble Machine learning Algorithms

Classifying transcripts as lncRNAs is difficult because there is no consensus
on what a lncRNA truly is.
Instead of using thresholds, or rules, for identifying lncRNAs, this tool uses
an ensemble stacking method of 8 different gradient boostling models to
predict lncRNAs.
Trained only on true, validated lncRNAs, this method has been tested on plant
transcriptomes with a high success rate.

See our publication at:
[Prediction of plant lncRNA by ensemble machine learning classifiers](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4665-2)

If you use CREMA, please cite:
Simopoulos CMA, Weretilnyk EA, Golding GB. "Prediction of plant lncRNA by ensemble machine learning classifiers." BMC Genomics (2018) doi:110.1186/s12864-018-4665-2

## Getting started

To use this tool, simply clone this repository on your machine by:
```
git clone https://github.com/gbgolding/crema.git
```

### Prerequisites

To use this tool you will need:

1. python3 
    - [biopython](http://biopython.org/)
    - [scikit-learn 0.20.0](http://scikit-learn.org)  
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

Firstly you must create the DIAMOND database from the [SwissProt protein database](http://www.uniprot.org/downloads):

```
diamond makedb --in uniprot_sprot.fasta -d swissprot.dmnd
```

Run DIAMOND:

```
diamond blastx -d swissprot.dmnd -q your_transcript_fasta_file.fa -o diamond_output.txt \\
-e 0.001 -k 5 --matrix BLOSUM62 --gapopen 11 --gapextend 1 --more-sensitive \\
-f 6 qseqid pident length qframe qstart qend sstart send evalue bitscore
```

Once you have identified your transcript features using CPAT and DIAMOND, you
can run the tool!

```
python3 bin/predict.py -f your_transcript_fasta_file.fa -c cpat_output.txt -d diamond_output.txt
```

Note: if script cannot find `logit_models.RData`, please run `predict.py` using its full file path. This is a known issue that we are working on solving.

## Output files

All output files are written to your working directory. Custom output directories to come...

The most helpful output file is `final_ensemble_predictions.csv`.
The CSV has outputs of both the features used in prediction as well as the lncRNA prediction score and final decision.

The columns describe:

1. gene name
2. length of transcript
3. ORF length
4. GC%
5. Fickett score (for more info see the [CPAT paper](https://academic.oup.com/nar/article/41/6/e74/2902455))
6. Hexamer score (for more info see the [CPAT paper](https://academic.oup.com/nar/article/41/6/e74/2902455))
7. % identity to a hit in the SwissProt database
8. Alignment length of hit in SwissProt database
9. Ratio of alignment length to transcript lenth
10. Ratio of alignment length to ORF length
11. Score of lncRNA prediction (you can use this to rank your predictions)
12. Final decision of prediction: 1 == lncRNA

The other files may be less useful to you, depending on what you're looking at.

`all_model_predictions.csv`: how each base model predicted the transcript (1 == lncRNA).   
`all_model_scores.csv`: the lncRNA prediction scores of each transcript for each base model.  
`ensemble_logreg_pred.csv`: the raw output of the final logistic regression stacking classifier.  

## Arguments

```
Required arguments:
    -f	input fasta file
    -c	output file from CPAT run
    -d	output file from Diamond blastx

Optional arguments:
    -s	minimum lncRNA prediction score (Default: 0.5)
```
