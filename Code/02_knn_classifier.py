#!/usr/bin/env python
# coding: utf-8

# In[81]:


import math
import pandas as pd
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import seaborn
import math
import operator
from scipy import stats
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import KNeighborsClassifier
import matplotlib.pyplot as plt
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score


# In[82]:


#reads in the data set
countsFile = "FRI data set - Gene counts.csv"
df= pd.read_csv(countsFile, dtype = {3: str})

#gets and drops other columns
metadata_cols = ['gene_name','gene_biotype', 'Chr', 'Start', 'End', 'Strand', 'Length']
metadata = df[metadata_cols]

#makes df with only needed columns
counts = df.drop(columns=metadata_cols)


# In[83]:


#subsets samples by acute and chronic
gene_count = len(counts)

#defines starting letter
start_letter_A = "A"

#creates list of columns that are acute
subset_cols_A = [col for col in df.columns if col.startswith(start_letter_A)]
subset_cols_A

#creates subsetted data frame for acute samples
metadata_cols_A = ['gene_id', 'gene_name']
acute_subset = df[metadata_cols_A + subset_cols_A]
acute_subset

#defines starting letter as c
start_letter_C = "C"

#creates list of columns that are chronic
subset_cols_C = [col for col in df.columns if col.startswith(start_letter_C)]
subset_cols_C

#creates subsetted data frame with chronic samples
metadata_cols_C = ['gene_id', 'gene_name']
chronic_subset = df[metadata_cols_C + subset_cols_C]


# In[84]:


#normalize data
metadata_1 = df.iloc[:, :8]
counts_1 = df.iloc[:, 8:]

#normalize (CPM)
cpm = counts_1.div(counts_1.sum(axis=0), axis = 1) * 1000000

#Filter data for top 10% of expressed genes


#log transformation
log2_cpm = np.log2(cpm + 1)

#makes acute groups
AP_cols = [c for c in log2_cpm.columns if c.startswith('AP')]
AE_cols = [c for c in log2_cpm.columns if c.startswith('AE')]

AP_data = log2_cpm[AP_cols]
AE_data = log2_cpm[AE_cols]

#log2foldchange
a_log2_fc = AE_data.mean(axis = 1) - AP_data.mean(axis = 1)

#t-test
a_t_stat, a_p_vals = stats.ttest_ind(AE_data, AP_data, axis = 1)

#makes table
results1 = metadata_1.copy()
results1['log2FoldChange'] = a_log2_fc
results1['p_value'] = a_p_vals

#filter out low counts
significant_genes1 = results1[(results1['p_value'] < 0.05) & (results1['log2FoldChange'].abs() > 1)]

#makes chronic groups
CP_cols = [c for c in log2_cpm.columns if c.startswith('CP')]
CE_cols = [c for c in log2_cpm.columns if c.startswith('CE')]

CP_data = log2_cpm[CP_cols]
CE_data = log2_cpm[CE_cols]

#log2foldchange
c_log2_fc = CE_data.mean(axis = 1) - CP_data.mean(axis = 1)

#t-test
c_t_stat, c_p_vals = stats.ttest_ind(CE_data, CP_data, axis = 1)

#makes table
results2 = metadata_1.copy()
results2['log2FoldChange'] = c_log2_fc
results2['p_value'] = c_p_vals

#filter out low counts
significant_genes2 = results2[(results2['p_value'] < 0.05) & (results2['log2FoldChange'].abs() > 1)]


# In[85]:


# Combine the significant gene indices from your two tests
significant_genes = list(set(significant_genes1.index) | set(significant_genes2.index))

# 2. Extract genes for AE and CE 
AE_subset = log2_cpm.loc[significant_genes, AE_cols].T  
CE_subset = log2_cpm.loc[significant_genes, CE_cols].T

# Assign labels for the classifier
AE_subset['label'] = 0  # Acute Experimental
CE_subset['label'] = 1  # Chronic Experimental

# Combine into a single dataframe
final_df = pd.concat([AE_subset, CE_subset])

# Define Features (X) and Target (y)
X = final_df.drop('label', axis=1)
y = final_df['label']


# In[86]:


# Updated split to 80% Training and 20% Testing
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, stratify=y, random_state=42)

# Scaling (Standardizing) based on the new 80% training set
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)


# In[87]:


# Classifier
knn = KNeighborsClassifier(n_neighbors=1)
knn.fit(X_train_scaled, y_train)

# Accuracy check
print(f"Model Accuracy: {knn.score(X_test_scaled, y_test)}")


# In[96]:


# 1. Define the range to go by 2s (start, stop, step)
k_range = range(1, 26, 2) 
test_accuracies = []

for k in k_range:
    knn = KNeighborsClassifier(n_neighbors=k)
    knn.fit(X_train_scaled, y_train)
    
    # Calculate test accuracy
    score = knn.score(X_test_scaled, y_test)
    test_accuracies.append(score * 100)

# 2. Create the plot
plt.figure(figsize=(10, 6))
plt.plot(k_range, test_accuracies, color='orange', marker='o', label='Test Accuracy')

# 3. Match the X-axis ticks to your K values
plt.xticks(k_range) 

plt.title('KNN: Testing Accuracy Acute vs Chronic')
plt.xlabel('K in KNN')
plt.ylabel('Accuracy (%)')
plt.grid(True, alpha=0.3)
plt.legend()

plt.show()


# In[ ]:




