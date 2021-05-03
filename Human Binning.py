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
# This file combines the human extracted data
# All human genes were combined into the same csv file before starting

data = pd.read_csv(r'combined_csv_file')
data


# In[3]:


data.drop_duplicates(subset ="Name",keep = 'first', inplace = True)
data


# In[4]:


data_GC3 = data.groupby(np.arange(len(data))//1000).mean()
data_GC3 = data_GC3.drop(columns=['TAA','TGA','TAG','TAAA','TAAT','TAAG','TAAC','TGAA','TGAT','TGAG','TGAC','TAGA','TAGT','TAGG','TAGC','Unk_Stop','Unk_four','sum'])
data_GC3


# In[5]:


data_stops = data.groupby(np.arange(len(data))//1000).sum()
data_stops = data_stops.drop(columns= ['Global_GC','Global_GC3','Unk_Stop','Unk_four','sum'])
data_stops


# In[6]:


all_data = pd.concat([data_GC3, data_stops], axis=1)
all_data


# In[7]:


data_stdev = data.groupby(np.arange(len(data))//1000).std()
data_stdev = data_stdev.drop(columns= ['Unk_Stop','Unk_four','sum'])
data_stdev


# In[8]:


data_all_stdev = pd.concat([all_data, data_stdev], axis=1)
data_all_stdev


# In[9]:


data_all_stdev.to_csv(r'exported.csv')


