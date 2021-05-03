#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import Bio as bio
import re
import os
import glob
from Bio import SeqIO
import pprint
from BioSQL import BioSeqDatabase
from operator import add


# In[2]:


i = 0
gb_file = r'file location.gbff'
mega_file = SeqIO.parse(open(gb_file,"r"), "genbank")
while i < 1:
	for record in mega_file:
    		i+=1
    		if i % 100 == 0:
        print(i)
print(i)


# # Code for extracting gene

# In[9]:


x = record.features[8]

y = str(x)
z = re.search('location:(.+?)]',y)
lo = z.group(1)
lo = lo[lo.find('['):]
if ' ' in lo:
    lo = lo.replace(' ','')
if '[' in lo:
    loc = lo.replace('[','')
mid = loc.find(':')
loc1 = int(loc[:mid])
loc2 = int(loc[mid+1:])
gene = str(record.seq[loc1:loc2])
print(record.features[1].type)
print(gene)
print(type(gene))



# # Function gene extraction returns gene

# In[11]:


def code_rev(data,j):
    x = data.features[j]
    y = str(x)
    z = re.search('location:(.+?)]',y)
    lo = z.group(1)
    loc = lo[lo.find('['):]
    if ' ' in loc:
        loc = loc.replace(' ','')
    if '[' in loc:
        loc = loc.replace('[','')
    if '{' in loc:
        loc = loc.replace('{','')
    mid = loc.find(':')
    loc1 = int(loc[:mid])
    loc2 = int(loc[mid+1:])
    gene = str(data.seq[loc1:loc2])
    if data.features[j].strand == -1:
        gene = gene[::-1]
        i = 0
        gene_rev = ''
        while i < len(gene):
            if gene[i] =='A':
                gene_rev +='T'
            elif gene[i] == 'T':
                gene_rev += 'A'
            elif gene[i] == 'G':
                gene_rev += 'C'
            elif gene[i] =='C':
                gene_rev += 'G'
            i+=1
        del(gene)
        gene = gene_rev
    return(gene)



# In[17]:


data = record



# In[19]:


def GC_3(gene):
    i=0
    stop = []
    GC3_NA = []
    while i+3 <= int((len(gene)/3)):
        if gene[(i*3+3)]== 'A': 
            GC3_NA.append('A')
            i +=1
        elif gene[(i*3+3)]== 'T': 
            GC3_NA.append('T')
            i+=1
        elif gene[(i*3+3)]== 'G': 
            GC3_NA.append('G')
            i+=1
        elif gene[(i*3+3)]== 'C': 
            GC3_NA.append('C')
            i+=1
        else:
            i+=1
    GC3 = 0
    GC3 = ((GC3_NA.count('G') + GC3_NA.count('C'))/(len(gene)/3))    
    return(GC3)


# In[20]:


def four_stop(data,j,codon_four):
    x = data.features[j]
    y = str(x)
    z = re.search('location:(.+?)]',y)
    lo = z.group(1)
    loc = lo[lo.find('['):]
    if ' ' in loc:
        loc = loc.replace(' ','')
    if '[' in loc:
        loc = loc.replace('[','')
    if '{' in loc:
        loc = loc.replace('{','')
    mid = loc.find(':')
    if data.features[j].strand ==1: 
        loc1 = int(loc[:mid])
        loc2 = int((loc[mid+1:]))+1
        gene = str(data.seq[loc1:loc2])
    elif data.features[j].strand == -1:
        loc1 = int(loc[:mid])-1
        loc2 = int((loc[mid+1:]))
        gene = str(data.seq[loc1:loc2])
        gene = gene[::-1]
        i = 0
        gene_rev = ''
        while i < len(gene):
            if gene[i] =='A':
                gene_rev +='T'
            elif gene[i] == 'T':
                gene_rev += 'A'
            elif gene[i] == 'G':
                gene_rev += 'C'
            elif gene[i] =='C':
                gene_rev += 'G'
            i+=1
        del(gene)
        gene = gene_rev
    stopper = gene[-4:]
    ar = np.array(codon_four)
    u = np.where(ar== stopper)
    try:
        p= [u[0]]
        q = int(p[0])
        codon_four[q][1] = codon_four[q][1]+1
    except:
        codon_four[12][1] = codon_four[12][1]+1
    return(codon_four)


# In[21]:


def file_runner(data):
    j = 0
    codon_data = [['TAA',0],['TAG',0],['TGA',0],['unk',0]]
    codon_unk =[]
    codon_GC_data = [['TAA'],['TAG'],['TGA']]
    codon_GC3_data = [['TAA'],['TAG'],['TGA']]
    global_gene = code_rev(data,0)
    global_GC = (global_gene.count('G')+global_gene.count('C'))/len(global_gene)
    global_GC3 = GC_3(global_gene)
    codon_four = [['TAAA',0],['TAAT',0],['TAAG',0],['TAAC',0],['TGAA',0],['TGAT',0],['TGAG',0],['TGAC',0],['TAGA',0],['TAGT',0],['TAGG',0],['TAGC',0],['unk',0]]
    while j < len(data.features):
        try:
            if data.features[j].type == 'CDS':
                gene = code_rev(data,j)
                codon = gene[-3:]
                GC_per = (gene.count('G')+gene.count('C'))/len(gene)
                GC_3per = GC_3(gene)
                if codon == 'TAA':
                    codon_data[0][1] = codon_data[0][1]+1
                    codon_GC_data[0].append(GC_per)
                    codon_GC3_data[0].append(GC_3per)
                    codon_four = four_stop(data,j,codon_four)
                    j+=1
                elif codon == 'TAG':
                    codon_data[1][1] = codon_data[1][1]+1
                    codon_GC_data[1].append(GC_per)
                    codon_GC3_data[1].append(GC_3per)
                    codon_four = four_stop(data,j,codon_four)
                    j+=1
                elif codon == 'TGA':
                    codon_data[2][1] = codon_data[2][1]+1
                    codon_GC_data[2].append(GC_per)
                    codon_GC3_data[2].append(GC_3per)
                    codon_four = four_stop(data,j,codon_four)
                    j+=1
                else:
                    codon_data[3][1] = codon_data[3][1]+1
                    codon_four = four_stop(data,j,codon_four)
                    j+=1
                    if len(codon_unk) ==0:
                        codon_unk.append(data.description)
                        codon_unk.append(codon)
                    else:
                        codon_unk.append(codon)
                    
            else:
                j+=1
        except:
            j+=1
    codon_df = pd.DataFrame(codon_GC_data, index = ['TAA','TAG','TGA'])
    codon_GC3df = pd.DataFrame(codon_GC3_data, index = ['TAA','TAG','TGA'])
    GC_list = (codon_df.mean(axis = 1))
    GC3_list = (codon_GC3df.mean(axis=1))
    Name_list = ['Name','TAA','TGA','TAG','Global_GC','Global_GC3','TAA_GC','TGA_GC','TAG_GC','TAA_GC3','TGA_GC3','TAG_GC3','TAAA','TAAT','TAAG','TAAC','TGAA','TGAT','TGAG','TGAC','TAGA','TAGT','TAGG','TAGC','Unk_Stop','Unk_four']
    file_data = [data.description,codon_data[0][1],codon_data[2][1],codon_data[1][1],(global_GC),(global_GC3),GC_list[0],GC_list[2],GC_list[1],GC3_list[0],GC3_list[2],GC3_list[1],codon_four[0][1],codon_four[1][1],codon_four[2][1],codon_four[3][1],codon_four[4][1],codon_four[5][1],codon_four[6][1],codon_four[7][1],codon_four[8][1],codon_four[9][1],codon_four[10][1],codon_four[11][1],codon_data[3][1],codon_four[12][1]]
    df_file_data = pd.DataFrame([file_data],columns = Name_list)
    df_file_data = df_file_data.set_index('Name')
    return([df_file_data,codon_unk])
    


# In[22]:


file_runner(record)


# In[23]:


Name_list = ['Name','TAA','TGA','TAG','Global_GC','Global_GC3','TAA_GC','TGA_GC','TAG_GC','TAA_GC3','TGA_GC3','TAG_GC3','TAAA','TAAT','TAAG','TAAC','TGAA','TGAT','TGAG','TGAC','TAGA','TAGT','TAGG','TAGC','Unk_Stop','Unk_four']


# In[24]:


archea_df = pd.DataFrame(columns = Name_list)
archea_df = archea_df.set_index('Name')
archea_df


# In[25]:


archea_df = pd.DataFrame(columns = Name_list)
archea_df = archea_df.set_index('Name')
holder_2 = file_runner(data)
full_archea_df = holder_2[0]
unknown_stops = pd.DataFrame(holder_2[1])
full_archea_df


# # Adjusted Runner each gene seperate
# 
# Only 10 repeats to check it is functioning

# In[251]:


i = 0
unk_stops = []
gb_file = r'C:\Users\Mole\Documents\Uni work\UNI WORK\Year 4\The project\Data\Human_RNA\human_1_rna.gbff'
mega_file = SeqIO.parse(open(gb_file,"r"), "genbank")
for record in mega_file:
    if i <10:
        holder = file_runner(record)
        archea_df = archea_df.append(holder[0])
        if holder[1] != []:
            unk_stops.append(holder[1])
        i+=1
    else:
        break
        
        


# In[214]:


archea_df


# In[26]:


archea_df = pd.DataFrame(columns = Name_list)
archea_df = archea_df.set_index('Name')
archea_df


# In[27]:
# Function for all repeats in a file

i = 0
unk_stops = []
gb_file = r'file_location.gbff'
mega_file = SeqIO.parse(open(gb_file,"r"), "genbank")
for record in mega_file:
    if i <10000000:
        holder = file_runner(record)
        archea_df = archea_df.append(holder[0])
        if holder[1] != []:
            unk_stops.append(holder[1])
        if i % 100 == 0:
            print(i)
        i+=1
    else:
        break
            
    


# In[277]:




mega_file


# In[28]:


archea_df


# In[112]:


unk_stops


# In[29]:


archea_df['sum'] = archea_df['TAA'] + archea_df['TGA'] + archea_df['TAG']
indexNames = archea_df[ archea_df['sum'] == 0 ].index
archea_df.drop(indexNames , inplace=True)

archea_df


# In[30]:


archea_df.to_csv(r'C:\Users\Mole\Documents\Uni work\UNI WORK\Year 4\The project\Extracted_csvs\human_8_gene_csv')



