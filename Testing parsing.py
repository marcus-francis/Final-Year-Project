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



# ## Code extracter reversed All combined into 1

# In[20]:


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


# ## Code Extracter 4 stop

# In[21]:


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


# ## GC 3 content

# In[22]:


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


# ## Intergene space explorerer

# In[23]:


def int_gene(j):
    x = data.features[j]
    y = str(x)
    z = re.search('location:(.+?)]',y)
    lo = z.group(1)
    if ' ' in lo:
        loc = lo.replace(' ','')
    if '[' in loc:
        loc = loc.replace('[','')
    if '{' in loc:
        loc = loc.replace('{','')
    mid = loc.find(':')
    loc1 = int(loc[:mid])
    loc2 = int(loc[mid+1:])
    if data.features[j].strand == -1:
        x = data.features[j-1]
        y = str(x)
        z = re.search('location:(.+?)]',y)
        lo = z.group(1)
        if ' ' in lo:
            lo = lo.replace(' ','')
            if '[' in lo:
                loc = lo.replace('[','')
        loce1 = int(loc[:mid])
        loce2 = int(loc[mid+1:])
        int_space = str(data.seq[loce2:loc1])
        int_space = int_space[::-1]
        i = 0
        int_rev = ''
        while i < len(int_space):
            if int_space[i] =='A':
                int_rev +='T'
            elif int_space[i] == 'T':
                int_rev += 'A'
            elif int_space[i] == 'G':
                int_rev += 'C'
            elif int_space[i] =='C':
                int_rev += 'G'
            i+=1
        del(int_space)
        int_space = int_rev
        return(int_space)
    else:
        x = data.features[j+1]
        y = str(x)
        z = re.search('location:(.+?)]',y)
        lo = z.group(1)
        if ' ' in lo:
            lo = lo.replace(' ','')
            if '[' in lo:
                loc = lo.replace('[','')
        loce1 = int(loc[:mid])
        loce2 = int(loc[mid+1:])
        int_space = str(data.seq[loc2:loce1])
        return(int_space)
        

    
    


# ## Extra stop codon finder

# In[24]:


def extra_stop(gene):
    i=0
    stop = []
    codon = []
    while i < int((len(gene)/3)):
        if gene[(i*3):(i*3+3)]== 'TAA': 
            stop.append(gene.find('TAA',i*3,i*3+3))
            codon.append('TAA')
        elif gene[(i*3):(i*3+3)]== 'TGA': 
            stop.append(gene.find('TGA',i*3,i*3+3))
            codon.append('TGA')
        elif gene[(i*3):(i*3+3)]== 'TAG': 
            stop.append(gene.find('TAG',i*3,i*3+3))
            codon.append('TAG')
        else:
            i+=1
    if codon != []:
        return(codon)
    else:
        pass

    
    


# ## Data Frame name

# In[25]:


Name_list = ['Name','TAA','TGA','TAG','Global_GC','Global_GC3','TAA_GC','TGA_GC','TAG_GC','TAA_GC3','TGA_GC3','TAG_GC3','TAAA','TAAT','TAAG','TAAC','TGAA','TGAT','TGAG','TGAC','TAGA','TAGT','TAGG','TAGC','Unk_Stop','Unk_four']


# ## Function of total file runner

# In[26]:


def race_runner(archea_df):
    files = (glob.glob(r'File location*.embl'))
    l = 0
    unk_stops = []
    while l < len(files):
        data = SeqIO.read(str(files[l]),'embl')
        holder = file_runner(data)
        archea_df = archea_df.append(holder[0])
        if holder[1] != []:
            unk_stops.append(holder[1])
        print(l,len(files))
        l+=1
    
    return([archea_df,unk_stops])


# ## Function file Runner

# In[27]:


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
    


# In[28]:


archea_df = pd.DataFrame(columns = Name_list)
archea_df = archea_df.set_index('Name')
archea_df


# In[1]:
## Full runner of files for extraction

archea_df = pd.DataFrame(columns = Name_list)
archea_df = archea_df.set_index('Name')
holder_2 = race_runner(archea_df)
full_archea_df = holder_2[0]
unknown_stops = pd.DataFrame(holder_2[1])
full_archea_df


# In[2]:


unknown_stops


# In[34]:


full_archea_df


# In[35]:


full_archea_df.to_csv(r'csv_file')


# In[36]:


unknown_stops.to_csv(r'csv_file')



