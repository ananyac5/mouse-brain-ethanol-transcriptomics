#!/usr/bin/env python
# coding: utf-8

# In[53]:


import math 
import pandas as pd
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import seaborn
from scipy import stats


# # Pre-Processing the Data

# In[54]:


#reads in the data set
countsFile = "FRI data set - Gene counts.csv"
df= pd.read_csv(countsFile, dtype = {3: str})

#gets and drops other columns
metadata_cols = ['gene_name','gene_biotype', 'Chr', 'Start', 'End', 'Strand', 'Length']
metadata = df[metadata_cols]

#makes df with only needed columns
counts = df.drop(columns=metadata_cols)


# In[55]:


#counts


# In[56]:


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


#chronic_subset


# # Pre-processing Data and finding DEG for Acute (treatment vs control)

# In[57]:


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

#prints statement
print(f"Found {len(significant_genes1)} significant genes.")
results1.to_csv('differential_expression_results.csv', index = False)

print(significant_genes1)


# # Pre-processing Data and finding DEG for Chronic (treatment vs control)

# In[58]:


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

#prints statement
print(f"Found {len(significant_genes2)} significant genes.")
significant_genes2.to_csv('c_differential_expression_results.csv', index = False)
print(significant_genes2)


# # Pre-processing Data and finding DEG for cortex

# In[27]:


#makes groups
ctx_control_cols = [c for c in log2_cpm.columns if c.startswith('CP_CTX') or c.startswith('AP_CTX')]
ctx_experiment_cols = [c for c in log2_cpm.columns if c.startswith('CE_CTX') or c.startswith('AE_CTX')]

ctx_control_data = log2_cpm[CTX_control_cols]
ctx_experiment_data = log2_cpm[CTX_experiment_cols]

#log2foldchange
ctx_log2_fc = CTX_experiment_data.mean(axis = 1) - CTX_control_data.mean(axis = 1)

#t-test
ctx_t_stat, ctx_p_vals = stats.ttest_ind(CTX_experiment_data, CTX_control_data, axis = 1)

#makes table
ctx_results = metadata_1.copy()
ctx_results['log2FoldChange'] = ctx_log2_fc
ctx_results['p_value'] = ctx_p_vals

#filter out low counts
ctx_significant_genes = ctx_results[(ctx_results['p_value'] < 0.05) & (ctx_results['log2FoldChange'].abs() > 1)]

#prints statement
print(f"Found {len(ctx_significant_genes)} significant genes.")
ctx_significant_genes.to_csv('ctx_differential_expression_results.csv', index = False)


# # Pre-processing Data and finding DEG for dorsal hippocampus

# In[29]:


#makes groups
dHPC_control_cols = [c for c in log2_cpm.columns if c.startswith('CP_dHPC') or c.startswith('AP_dHPC')]
dHPC_experiment_cols = [c for c in log2_cpm.columns if c.startswith('CE_dHPC') or c.startswith('AE_dHPC')]

dHPC_control_data = log2_cpm[dHPC_control_cols]
dHPC_experiment_data = log2_cpm[dHPC_experiment_cols]

#log2foldchange
dHPC_log2_fc = dHPC_experiment_data.mean(axis = 1) - dHPC_control_data.mean(axis = 1)

#t-test
dHPC_t_stat, dHPC_p_vals = stats.ttest_ind(dHPC_experiment_data, dHPC_control_data, axis = 1)

#makes table
dHPC_results = metadata_1.copy()
dHPC_results['log2FoldChange'] = dHPC_log2_fc
dHPC_results['p_value'] = dHPC_p_vals

#filter out low counts
dHPC_significant_genes = dHPC_results[(dHPC_results['p_value'] < 0.05) & (dHPC_results['log2FoldChange'].abs() > 1)]

#prints statement
print(f"Found {len(dHPC_significant_genes)} significant genes.")
dHPC_significant_genes.to_csv('dHPC_differential_expression_results.csv', index = False)


# # Pre-processing Data and finding DEG for ventral hippocampus

# In[36]:


#makes groups
vHPC_control_cols = [c for c in log2_cpm.columns if c.startswith('CP_vHPC') or c.startswith('AP_vHPC')]
vHPC_experiment_cols = [c for c in log2_cpm.columns if c.startswith('CE_vHPC') or c.startswith('AE_vHPC')]

vHPC_control_data = log2_cpm[vHPC_control_cols]
vHPC_experiment_data = log2_cpm[vHPC_experiment_cols]

#log2foldchange
vHPC_log2_fc = vHPC_experiment_data.mean(axis = 1) - vHPC_control_data.mean(axis = 1)

#t-test
vHPC_t_stat, vHPC_p_vals = stats.ttest_ind(vHPC_experiment_data, vHPC_control_data, axis = 1)

#makes table
vHPC_results = metadata_1.copy()
vHPC_results['log2FoldChange'] = vHPC_log2_fc
vHPC_results['p_value'] = vHPC_p_vals

#filter out low counts
vHPC_significant_genes = vHPC_results[(vHPC_results['p_value'] < 0.05) & (vHPC_results['log2FoldChange'].abs() > 1)]

#prints statement
print(f"Found {len(vHPC_significant_genes)} significant genes.")
vHPC_significant_genes.to_csv('vHPC_differential_expression_results.csv', index = False)


# # Pre-processing Data and finding DEG for dorsal striatum

# In[38]:


#makes groups
DS_control_cols = [c for c in log2_cpm.columns if c.startswith('CP_DS') or c.startswith('AP_DS')]
DS_experiment_cols = [c for c in log2_cpm.columns if c.startswith('CE_DS') or c.startswith('AE_DS')]

DS_control_data = log2_cpm[DS_control_cols]
DS_experiment_data = log2_cpm[DS_experiment_cols]

#log2foldchange
DS_log2_fc = DS_experiment_data.mean(axis = 1) - DS_control_data.mean(axis = 1)

#t-test
DS_t_stat, DS_p_vals = stats.ttest_ind(DS_experiment_data, DS_control_data, axis = 1)

#makes table
DS_results = metadata_1.copy()
DS_results['log2FoldChange'] = DS_log2_fc
DS_results['p_value'] = DS_p_vals

#filter out low counts
DS_significant_genes = DS_results[(DS_results['p_value'] < 0.05) & (DS_results['log2FoldChange'].abs() > 1)]

#prints statement
print(f"Found {len(DS_significant_genes)} significant genes.")
DS_significant_genes.to_csv('DS_differential_expression_results.csv', index = False)


# # Pre-processing Data and finding DEG for ventral striatum

# In[41]:


#makes groups
VS_control_cols = [c for c in log2_cpm.columns if c.startswith('CP_VS') or c.startswith('AP_VS')]
VS_experiment_cols = [c for c in log2_cpm.columns if c.startswith('CE_VS') or c.startswith('AE_VS')]

VS_control_data = log2_cpm[VS_control_cols]
VS_experiment_data = log2_cpm[VS_experiment_cols]

#log2foldchange
VS_log2_fc = VS_experiment_data.mean(axis = 1) - VS_control_data.mean(axis = 1)

#t-test
VS_t_stat, VS_p_vals = stats.ttest_ind(VS_experiment_data, VS_control_data, axis = 1)

#makes table
VS_results = metadata_1.copy()
VS_results['log2FoldChange'] = VS_log2_fc
VS_results['p_value'] = VS_p_vals

#filter out low counts
VS_significant_genes = VS_results[(VS_results['p_value'] < 0.05) & (VS_results['log2FoldChange'].abs() > 1)]

#prints statement
print(f"Found {len(VS_significant_genes)} significant genes.")
VS_significant_genes.to_csv('VS_differential_expression_results.csv', index = False)


# # Heat Map Visualization for Accute vs Chronic

# In[61]:


import matplotlib.pyplot as plt
import numpy as np

def volcano_plot(results, title):
    plt.figure(figsize=(8, 6))
    
    # define colors
    colors = ['red' if (p < 0.05 and lfc > 1) else 'blue' if (p < 0.05 and lfc < -1) else 'grey'
              for lfc, p in zip(results['log2FoldChange'], results['p_value'])]
    
    plt.scatter(results['log2FoldChange'], -np.log10(results['p_value']), 
                c=colors, alpha=0.5, s=10)
    
    # threshold lines
    plt.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=0.8)
    plt.axvline(1, color='black', linestyle='--', linewidth=0.8)
    plt.axvline(-1, color='black', linestyle='--', linewidth=0.8)
    
    plt.xlabel('log2 Fold Change')
    plt.ylabel('-log10(p-value)')
    plt.title(title)
    plt.show()

# Acute
#volcano_plot(results1, 'Volcano Plot — Acute Ethanol vs Placebo')

# Chronic
#volcano_plot(results2, 'Volcano Plot — Chronic Ethanol vs Placebo')

# CTX
#volcano_plot(ctx_results, 'Volcano Plot — Ethanol in Cortex vs Placebo')

# dHPC
#volcano_plot(dHPC_results, 'Volcano Plot — Ethanol in Dorsal Hippocampus vs Placebo')

# vHPC
#volcano_plot(vHPC_results, 'Volcano Plot — Ethanol in Ventral Hippocampus vs Placebo')

# DS
volcano_plot(DS_results, 'Volcano Plot — Ethanol in Dorsal Striatum vs Placebo')

# VS
#volcano_plot(VS_results, 'Volcano Plot — Ethanol in Ventral Striatum vs Placebo')


# In[14]:


significant_indices = list(set(significant_genes1.index) | set(significant_genes2.index))
DEG = log2_cpm.loc[significant_indices]
heat_map = seaborn.clustermap(
    DEG, 
    cmap="viridis", 
    z_score=0,            # Standardizes data across the rows (genes)
    figsize=(10, 8),      # Adjusts the size of the plot
    yticklabels=False     # Hides gene names on the y-axis (useful if you have hundreds of genes)
)

# Add a title and display
heat_map.figure.suptitle('Expression of Significant DEGs (Acute & Chronic)', y=1.05)
plt.show()


# In[ ]:





# In[ ]:





# In[ ]:




