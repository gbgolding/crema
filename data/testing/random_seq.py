from Bio import SeqIO
from random import sample
import csv

with open("/home/caitlin/lncRNApred/known_lncRNAs/gb/testing_models/athal/TAIR10_cdna_20101214_updated_CLEANED.fa") as f:
    seqs = list(SeqIO.parse(f, 'fasta'))

random = sample(seqs, 100)

random_names = []
for record in random:
    random_names.append(record.id)

random_athal = open('athaliana_random100.fa', 'w')

SeqIO.write(random, random_athal, 'fasta')

with open('athaliana_random100_names.txt', 'w') as f:
    a = csv.writer(f, delimiter=',')
    a.writerows(random_names)
