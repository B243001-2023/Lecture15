# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 10:57:41 2023

@author: s2600370
"""

#%%
#Exercise of dict
enzymes = { 
'EcoRI' : 'GAATTC',
'AvaII' : 'GGACC',
'BisI'  : 'GCNGC' 
}

enzymes['EcoRI']


#%%
dna = 'ATCGTATCGATGTACGCTGA'

dna_dict = dict()
for i in set(dna):
    dna_dict[i] = dna.count(i)


#%%
c = dict(i=1 ,k=2)


#%%
enzymes['EcoRI'] = 'GAATTC'
enzymes['AvaII'] = 'GGACC'
enzymes['BisI'] = 'GCNGC'
enzymes['RanDomI'] = 'GCNYC'
enzymes['RanDomI'] = 'NewYorkCITY'
enzymes['RanDomI'] = enzymes['RanDomI'] + ' ABCDE,34'

enzymes2 = { 'EcoRI' : 'GAATTC','AvaII' : 'GGACC','BisI' : 'GCNGC' }

#%%
dna = "AATGATGAACGAC"
dinucleotides=[]
for first in 'ATGC' :
  for second in 'ATGC' :
    dinucleotides.append(str(first)+str(second))

dinucleotides
all_counts = dict()

for dinucleotide in dinucleotides: 
    count = dna.upper().count(dinucleotide) 
    print("count is " + str(count) + " for " + dinucleotide) 
    all_counts[dinucleotide] = count 
    
#removing the 0
all_counts2 = {key: value for key, value in all_counts.items() if value != 0}

#%% exercise
#q1
qs = [
"What's your name?",
"How old are you?",
"What is your favourite colour?",
"Do you like Python?",
"The world is flat: True or False?"
]

ans = dict()

for i in qs:
    ans[i] = input(i)
    print(ans[i])

#q2
gencode = {
'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

dna = input('please input your seq')
protein = ['','','']
for format in [0,1,2]:
    for i in range(format, len(dna), 3):
        code = dna[i:i+3].upper()
        if 'N' in code:
            protein[format] += '?'
        elif len(code)==3:
            protein[format] += gencode[code]
            code = ''
print(protein)

rev_dna = {'A':'T', 'G':'C', 'C':'G', 'T':'A'}
gencode2 = dict()
rev_key = ''
for i in gencode.keys():
    rev_key = rev_dna[i[0]] + rev_dna[i[1]] + rev_dna[i[2]]
    gencode2[rev_key] = gencode[i]
   
rev_protein = ['','','']
for format in [0,1,2]:
    for i in range(format, len(dna), 3):
        code = dna[i:i+3].upper()
        if 'N' in code:
            rev_protein[format] += '?'
        elif len(code)==3:
            rev_protein[format] += gencode2[code]
            code = ''
print(rev_protein)
  
    
  
    
  
    