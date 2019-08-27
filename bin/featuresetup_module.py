#!/usr/bin/python

from Bio import SeqIO
from Bio.SeqUtils import GC
from sklearn import svm
import numpy as np
import csv
import operator
import os 

## use previously subsetting dataset for ease


def transcript_info_dict(fasta_file, cpat_file, blast_file): #file should be in quotes
    name_list = []
    transcript_dict = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        name = record.id
        name = name.lower()
        transcript_dict[name] = {}
        transcript_dict[name]["length"] = len(record.seq)   
        transcript_dict[name]["GC"] = GC(record.seq)    
        name_list.append(name.lower())
    with open(cpat_file, 'r') as tabular:
        cpat_reader = csv.reader(tabular, delimiter='\t')
        next(cpat_reader, None) # skip header
        for row in cpat_reader:
            name = row[0]
            name = name.lower() 
            ORF = float(row[2]) 
            score = float(row[5])
            fickett = float(row[3])
            hexamer = float(row[4])
            transcript_dict[name]["ORF"] = ORF  
            transcript_dict[name]["fickett"] = fickett 
            transcript_dict[name]["hexamer"] = hexamer  
            #transcript_dict[name]["coding_prob"] = score # we are going to take out coding_prob, and just use hexamer and fickett scores...the parameters of the CPAT model
    ## Appends scores to "maximum" until new gene name is found
    print(list(transcript_dict)[0]) 
    if os.stat(blast_file).st_size == 0:
        for gene in name_list: 
            transcript_dict[gene]["identity"] = float(0)
            transcript_dict[gene]["align_length"] = float(0)
            transcript_dict[gene]["align_perc_len"] = 0
            transcript_dict[gene]["align_perc_ORF"] = 0         
    else:    
        with open(blast_file, "r") as f:
            tab_reader = csv.reader(f, delimiter='\t')    
            line_1 = next(tab_reader)
            first = line_1[0]
            score = [float(line_1[9])]
            with_len = [[first, float(line_1[1]), float(line_1[2]), float(line_1[3]), float(line_1[9])]] # name identity length frame score 
            for row in tab_reader:
                if row[0] == first:
                    score.append(float(row[9]))
                    with_len.append([row[0], float(row[1]), float(row[2]), float(row[3]), float(row[9])])
                else:
                    max_value = max(score)
                    max_index = score.index(max_value) # should techincally all be '0', but just making sure
                    max_len_ident = with_len[max_index]
    #               max_score_len_iden.append(max_len_ident)
                    if max_len_ident[3] > 0: # only considers a hit if the first blast match is in positive frame
                        transcript_dict[first.lower()]["identity"] = float(max_len_ident[1])
                        transcript_dict[first.lower()]["align_length"] = float(max_len_ident[2])
                        transcript_dict[first.lower()]["align_perc_len"] = float(transcript_dict[first.lower()]["align_length"])/float(transcript_dict[first.lower()]["length"])
                        if float(transcript_dict[first.lower()]["ORF"]) > 0:
                            transcript_dict[first.lower()]["align_perc_ORF"] = float(transcript_dict[first.lower()]["align_length"])/float(transcript_dict[first.lower()]["ORF"])   
                        else:
                            transcript_dict[first.lower()]["align_perc_ORF"] = float(0)
                    score = [float(row[9])]
                    first = row[0]
                    with_len = [[first, float(row[1]), float(row[2]), float(row[3]), float(row[9])]] # set up original with_len list
        ## add 0 identity and 0 length to genes without blast hits!!!!
        for gene in name_list:
            if gene in transcript_dict and not "identity" in transcript_dict[gene]:
                transcript_dict[gene]["identity"] = float(0)
                transcript_dict[gene]["align_length"] = float(0)
                transcript_dict[gene]["align_perc_len"] = 0
                transcript_dict[gene]["align_perc_ORF"] = 0         
    #print(transcript_dict)
    # convert dictionary to np.array
    transcript_info_array = np.array([[transcript_dict[gene][feature] for feature in sorted(transcript_dict[gene])] for gene in name_list], dtype=float)
    #transcript_info_array.astype(float)    
    return(transcript_info_array, transcript_dict, name_list) 


def transcript_info(fasta_file, cpat_file, blast_file): #file should be in quotes
    name_list = []
    transcript_dict = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        name = record.id
        name = name.lower()
        transcript_dict[name] = {}
        transcript_dict[name]["length"] = len(record.seq)   
        transcript_dict[name]["GC"] = GC(record.seq)        
        name_list.append(name.lower())
    with open(cpat_file, "r") as tabular:
        cpat_reader = csv.reader(tabular, delimiter=("\t"))
        next(cpat_reader, None) # skip header
        for row in cpat_reader:
            name = row[0]
            name = name.lower() 
            ORF = float(row[2]) 
            score = float(row[5])
            fickett = float(row[3])
            hexamer = float(row[4])
            transcript_dict[name]["ORF"] = ORF  
            transcript_dict[name]["fickett"] = fickett 
            transcript_dict[name]["hexamer"] = hexamer  
            #transcript_dict[name]["coding_prob"] = score # we are going to take out coding_prob, and just use hexamer and fickett scores...the parameters of the CPAT model
    ## Appends scores to "maximum" until new gene name is found
    with open(blast_file, "r") as f:
        tab_reader = csv.reader(f, delimiter=("\t"))    
        line_1 = next(tab_reader)
        first = line_1[0]
        score = [float(line_1[9])]
        with_len = [[first, float(line_1[1]), float(line_1[2]), float(line_1[3]), float(line_1[9])]] # name identity length frame score 
        for row in tab_reader:
            if row[0] == first:
                score.append(float(row[9]))
                with_len.append([row[0], float(row[1]), float(row[2]), float(row[3]), float(row[9])])
            else:
                max_value = max(score)
                max_index = score.index(max_value) # should techincally all be '0', but just making sure
                max_len_ident = with_len[max_index]
#               max_score_len_iden.append(max_len_ident)
                if max_len_ident[3] > 0: # only considers a hit if the first blast match is in positive frame
                    transcript_dict[first.lower()]["identity"] = float(max_len_ident[1])
                    transcript_dict[first.lower()]["align_length"] = float(max_len_ident[2])
                    transcript_dict[first.lower()]["align_perc_len"] = float(transcript_dict[first.lower()]["align_length"])/float(transcript_dict[first.lower()]["length"])
                    if float(transcript_dict[first.lower()]["ORF"]) > 0:
                        transcript_dict[first.lower()]["align_perc_ORF"] = float(transcript_dict[first.lower()]["align_length"])/float(transcript_dict[first.lower()]["ORF"])   
                    else:
                        transcript_dict[first.lower()]["align_perc_ORF"] = 0
                score = [float(row[9])]
                first = row[0]
                with_len = [[first, float(row[1]), float(row[2]), float(row[3]), float(row[9])]] # set up original with_len list
    ## add 0 identity and 0 length to genes without blast hits!!!!
    for gene in name_list:
        if gene in transcript_dict and not "identity" in transcript_dict[gene]:
            transcript_dict[gene]["identity"] = float(0)
            transcript_dict[gene]["align_length"] = float(0)
            transcript_dict[gene]["align_perc_len"] = 0
            transcript_dict[gene]["align_perc_ORF"] = 0 
    # convert dictionary to np.array
    transcript_info_array = np.array([[transcript_dict[gene][feature] for feature in sorted(transcript_dict[gene])] for gene in name_list],dtype=float)
    return(transcript_info_array, name_list) 


### need to TRAIN not on COOLAIR and COLDAIR when testing for them!!!!!
def trans_info_dict_cc(fasta_file, cpat_file, blast_file): #file should be in quotes
    name_list = []
    transcript_dict = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        name = record.id
        name = name.lower()
        transcript_dict[name] = {}
        transcript_dict[name]["length"] = len(record.seq)   
        transcript_dict[name]["GC"] = GC(record.seq)    
        name_list.append(name.lower())
    with open(cpat_file, "r") as tabular:
        cpat_reader = csv.reader(tabular, delimiter=("\t"))
        next(cpat_reader, None) # skip header
        for row in cpat_reader:
            name = row[0]
            name = name.lower() 
            ORF = float(row[2]) 
            score = float(row[5])
            fickett = float(row[3])
            hexamer = float(row[4])
            transcript_dict[name]["ORF"] = ORF  
            transcript_dict[name]["fickett"] = fickett 
            transcript_dict[name]["hexamer"] = hexamer  
            #transcript_dict[name]["coding_prob"] = score # we are going to take out coding_prob, and just use hexamer and fickett scores...the parameters of the CPAT model
    ## Appends scores to "maximum" until new gene name is found
    with open(blast_file, "r") as f:
        tab_reader = csv.reader(f, delimiter=("\t"))    
        line_1 = next(tab_reader)
        first = line_1[0]
        score = [float(line_1[9])]
        with_len = [[first, float(line_1[1]), float(line_1[2]), float(line_1[3]), float(line_1[9])]] # name identity length frame score 
        for row in tab_reader:
            if row[0] == first:
                score.append(float(row[9]))
                with_len.append([row[0], float(row[1]), float(row[2]), float(row[3]), float(row[9])])
            else:
                max_value = max(score)
                max_index = score.index(max_value) # should techincally all be '0', but just making sure
                max_len_ident = with_len[max_index]
#               max_score_len_iden.append(max_len_ident)
                if max_len_ident[3] > 0: # only considers a hit if the first blast match is in positive frame
                    transcript_dict[first.lower()]["identity"] = float(max_len_ident[1])
                    transcript_dict[first.lower()]["align_length"] = float(max_len_ident[2])
                    transcript_dict[first.lower()]["align_perc_len"] = float(transcript_dict[first.lower()]["align_length"])/float(transcript_dict[first.lower()]["length"])
                    if float(transcript_dict[first.lower()]["ORF"]) > 0:
                        transcript_dict[first.lower()]["align_perc_ORF"] = float(transcript_dict[first.lower()]["align_length"])/float(transcript_dict[first.lower()]["ORF"])   
                    else:
                        transcript_dict[first.lower()]["align_perc_ORF"] = float(0)
                score = [float(row[9])]
                first = row[0]
                with_len = [[first, float(row[1]), float(row[2]), float(row[3]), float(row[9])]] # set up original with_len list
    ## add 0 identity and 0 length to genes without blast hits!!!!
    for gene in name_list:
        if gene in transcript_dict and not "identity" in transcript_dict[gene]:
            transcript_dict[gene]["identity"] = float(0)
            transcript_dict[gene]["align_length"] = float(0)
            transcript_dict[gene]["align_perc_len"] = 0
            transcript_dict[gene]["align_perc_ORF"] = 0         
    # remove coolair and coldair from transcript_dict AND name_list
    to_remove = ['coolair_arabidopsisthaliana_1', 'coldair_arabidopsisthaliana_1']
    for k in to_remove:
        transcript_dict.pop(k, None)
    name_list = [v for v in name_list if v not in to_remove]    
    # convert dictionary to np.array
    transcript_info_array = np.array([[transcript_dict[gene][feature] for feature in sorted(transcript_dict[gene])] for gene in name_list], dtype=float)
    #transcript_info_array.astype(float)    
    return(transcript_info_array, transcript_dict, name_list) 




## testing function

#trans_array, trans_names = transcript_info("subset.fa", "subset_cpat.txt", "subset.fa.tab") 
# both variables saved!


