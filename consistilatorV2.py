##########################################################
## Consistilator:  consistilatorV2.py                   ##
##  ______     ______     __  __                        ##
## /\  __ \   /\  ___\   /\ \/\ \                       ##
## \ \  __ \  \ \___  \  \ \ \_\ \                      ##
##  \ \_\ \_\  \/\_____\  \ \_____\                     ##
##   \/_/\/_/   \/_____/   \/_____/                     ##
## @Developed by: Plaisier Lab                          ##
##   (https://plaisierlab.engineering.asu.edu/)         ##
##   Arizona State University                           ##
##   242 ISTB1, 550 E Orange St                         ##
##   Tempe, AZ  85281                                   ##
## @Author:  Chris Plaisier                             ##
## @License:  GNU GPLv3                                 ##
##                                                      ##
## If this program is used in your analysis please      ##
## mention who built it. Thanks. :-)                    ##
##########################################################

from subprocess import *
from multiprocessing import Pool, cpu_count, Manager
from scipy.stats import pearsonr
import numpy as np
import pandas as pd
import statsmodels.api as sm
import json
import argparse

parser = argparse.ArgumentParser(description='use for testing argparse.py')
parser.add_argument('--tumor', help='tumor type', type = str)
args = parser.parse_args()


# Get a correlation p-value from R
def correlation(a1, a2):
    """
    Calculate the correlation coefficient and p-value between two variables.
    Input: Two arrays of float or integers.
    Returns: Corrleation coefficient and p-value.
    """
    """
    # Fire up R
    rProc = Popen('R --no-save --slave', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    runMe = []
    # Make the data into an R matrix
    runMe.append('c1 = cor.test(c('+','.join([str(i) for i in a1])+'),c('+','.join([str(i) for i in a2])+'))')
    runMe.append('c1$estimate')
    runMe.append('c1$p.value')
    runMe = '\n'.join(runMe)+'\n'
    out = rProc.communicate(runMe)
    # Process output
    splitUp = out[0].strip().split('\n')
    rho = float(splitUp[1])
    pValue = float((splitUp[2].split(' '))[1])
    """
    r1 = pearsonr(a1,a2)
    return [r1[0], r1[1]]

# Send model for comparison versus
def regression_R(response, predictors):
    """
    Fit a linear regression model of all terms.
    Input: A response variable and a specified number of predictor variables.
        - response = vector of floats
        - predictors = dictionary of vectors of floats hashed on variable names
    Returns: Overall model significance and term significance.
    """
    # Fire up R
    rProc = Popen('R --no-save --slave', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    runMe = []
    # Make the data into an R matrix
    tmp = []
    names = {}
    revNames = {}
    for i in predictors:
        names[i] = i.replace('-','_')
        revNames[i.replace('-','_')] = i
        tmp.append(names[i]+' = c('+','.join([str(i) for i in predictors[i]])+')')
    runMe.append('d1 = data.frame(response = c('+','.join([str(i) for i in response])+')'+','+','.join(tmp)+')')
    runMe.append('slm1 = summary(lm(response ~ .,data=d1))')
    runMe.append('slm1$coefficients')
    runMe.append('slm1$adj.r.squared')
    runMe = '\n'.join(runMe)+'\n'
    out = rProc.communicate(runMe)
    splitUp = out[0].strip().split('\n')
    adj_r_squared = splitUp.pop(-1).split(' ')[1]
    splitUp.pop(0)
    splitUp.pop(0)
    splitUp = [[j for j in i.split(' ') if j] for i in splitUp]
    tmp = dict(zip([revNames[i[0]] for i in splitUp],[{'estimate':i[1],'stdErr':i[2],'t':i[3],'p_value':i[4]} for i in splitUp]))
    tmp['adj_r_squared'] = adj_r_squared
    return tmp

# Send model for comparison versus
def regression(response, predictors):
    """
    Fit a linear regression model of all terms.
    Input: A response variable and a specified number of predictor variables.
        - response = a Pandas dataframe of the response
        - predictors = a Pandas dataframe of predictors
    Returns: Overall model significance and term significance.
    """
    # Run linear regression in Python
    tmp = {}
    model = sm.OLS(response, predictors).fit()
    for p1 in list(predictors.columns.values):
        tmp[p1] = {'estimate':model.params[p1],'stdErr':model.bse[p1],'t':model.tvalues[p1],'p_value':model.pvalues[p1]}
    tmp['adj_r_squared'] = model.rsquared_adj
    return tmp

def dfMe(obj):
    if isinstance(obj,pd.DataFrame):
        return obj.transpose()
    elif isinstance(obj,pd.Series):
        return obj.to_frame()

### Load sample subsets ###
#subsets = { 'all':range(0,12) }

# read in expression data
expression = {}
for gse in ['all']:
    expression[gse] = pd.read_csv('exprs_all/'+args.tumor+'_RNAseq.csv', header=0, index_col=0).dropna(axis=1) #pd.read_csv('tfExp.csv', header=0, index_col=0).transpose()

# Structure of Studer datasets
subsets = { 'all':list(expression['all'].columns) }

### Load sub-networks ###
# Load up FANMOD subgraph enumeration results
inFile = open('FanmodOutput/fanmod_motifs_'+args.tumor+'.txt','r')
# Get rid of headers
while 1:
    line = inFile.readline()
    if line.strip()=='Result overview:':
        inFile.readline()
        inFile.readline()
        inFile.readline()
        inFile.readline()
        break
subnetworks = {}
filter_pv = 0.05
filter_z = 2
while 1:
    line = inFile.readline()
    if not line:
        break
    line2 = inFile.readline().strip() # Get rid of second adjacency matrix line
    line3 = inFile.readline().strip() # Get rid of third adjacency matrix line
    #line4 = inFile.readline().strip() # Get rid of fourth adjacency matrix line
    #print line4
    inFile.readline() # Get rid of blank line
    splitUp = [i for i in line.strip().split(' ') if i]
    id = splitUp[1]+line2+line3 #+line4
    if float(splitUp[6]) <= filter_pv and float(splitUp[5]) >= filter_z:
        subnetworks[id] = splitUp[0]
inFile.close()

# Take out extrodinarily large subnetworks
#subnetworks = [i for i in subnetworks if not i in [74, 2184]]

# Create dictionary to convert
id2gene = {}
inFile = open('PanOutput/mdraw_'+args.tumor+'_mdraw_LIST.txt','r')
while 1:
    line = inFile.readline()
    if not line:
        break
    splitUp = line.strip().split(' ')
    id2gene[splitUp[1]] = splitUp[0]
inFile.close()

# Load up network motifs
networkMotifs = {}
inFile = open('FanmodOutput/fanmod_motifs_'+args.tumor+'.txt.dump','r')
inFile.readline() # Get rid of header
inFile.readline() # Get rid of header
while 1:
    line = inFile.readline()
    if not line:
        break
    splitUp = line.strip().split(',')
    #numMotif = int(splitUp.pop(0).replace('2','1'), 2)
    motif = splitUp.pop(0)
    if motif in subnetworks:
        if not motif in networkMotifs:
            networkMotifs[motif]  = []
        networkMotifs[motif].append([id2gene[i] for i in splitUp])
inFile.close()

# Read in Biotapesty file
outEdges = {}
inEdges = {}
inFile = open('PanOutput/biotapestry_CHIR_curve'+args.tumor+'.csv','r')
counts = 0
while 1:
    line = inFile.readline()
    if not line:
        break
    splitUp = line.strip().split(',')
    if splitUp[0]=='"# Standard Interactions"':
        inFile.readline() # Remove header
        break
while 1:
    line = inFile.readline()
    if not line:
        break
    splitUp = line.strip().split(',')
    node1 = splitUp[3].strip('"')
    node2 = splitUp[5].strip('"')
    if not node1 in outEdges:
        outEdges[node1] = []
    if not node2 in outEdges[node1]:
        outEdges[node1].append(node2)
    if not node2 in inEdges:
        inEdges[node2] = []
    if not node1 in inEdges[node2]:
        inEdges[node2].append(node1)
inFile.close()

# To translate gene ids later
symbol2entrez = {}
entrez2symbol = {}
with open('id_conversion/gene2entrezId.csv','r') as inFile:
    while 1:
        inLine = inFile.readline()
        if not inLine:
            break
        split = inLine.strip().split(',')
        entrez2symbol[split[1]] = split[0]
        symbol2entrez[split[0]] = split[1]


# To test the survival function
#print regression(expression['Cluster-10'], {'SP8':expression['SP8'], 'SIM2':expression['SIM2'],'NKX2-1':expression['NKX2-1']})

def consistentInstance(instance, inEdges_sm, expression_sm, writeMe, cur, symbol2entrez):
    print '    '+' '.join(instance)+' - '+cur['subset']+' - '+str(cur['subnetwork'])
    entry = {}
    for node in instance:
        # A. Determine inputs for each node
        inputs = []
        if node in inEdges_sm:
            inputs = [i for i in list(set(inEdges_sm[node]).intersection(instance)) if not i==node]
        # B. Test model consistency for inputs to each node
        if len(inputs)>0 and len([i for i in inputs+[node] if not int(symbol2entrez[i]) in list(expression_sm[cur['subset']].index)])==0:
            res1 = regression(dfMe(expression_sm[cur['subset']].loc[int(symbol2entrez[node])]),dfMe(expression_sm[cur['subset']].loc[[int(symbol2entrez[i]) for i in inputs]]))
            consistent = 'Inconsistent'
            if len([i for i in inputs if float(res1[int(symbol2entrez[i])]['p_value'])<=0.05])==len(inputs):
                consistent = res1['adj_r_squared']
            entry[node] = {'inputs':inputs,'consistency':consistent}
        else:
            entry[node] = {'inputs':inputs,'consistency':'NA'}
    writeMe.append(cur['subset']+','+cur['subnetwork_num']+','+str(cur['subnetwork'])+','+';'.join(instance)+','+','.join([i+','+';'.join(entry[i]['inputs'])+','+str(entry[i]['consistency']) for i in instance]))

# Make shared memory objects
#cpus = cpu_count()
#mgr = Manager()
#inEdges_sm = mgr.dict(inEdges)
inEdges_sm = inEdges
#expression_sm = mgr.dict(expression)
expression_sm = expression

### For each subset of samples including all ###
# Goal is to write out a file with this header:
# Data Subset,Netowrk Motif,Instance,Node1,Node1.Inputs,Node1.Consistency,Node2,Node2.Inputs,Node2.Consistency,Node3,Node3.Inputs,Node3.Consistency,Node4,Node4.Inputs,Node4.Consistency
for subset in ['all']:
    writeMe = [] #mgr.list()
    writeMe.append('Data Subset,Netowrk Motif,Adjacency Matrix,Instance,Node1,Node1.Inputs,Node1.Consistency,Node2,Node2.Inputs,Node2.Consistency,Node3,Node3.Inputs,Node3.Consistency')
    print 'Working on '+subset+' data subset.'
    ### For each network motif enriched in network ###
    for subnetwork in subnetworks:
        print '  Working on '+str(subnetwork)+' subnetwork.'
        ### For each subnetwork instance ###
        for instance in networkMotifs[subnetwork]:
            cur = {'subset':subset, 'samples': subsets[subset], 'subnetwork':subnetwork, 'subnetwork_num':subnetworks[subnetwork]}
            consistentInstance(instance, inEdges_sm, expression_sm, writeMe, cur, symbol2entrez)
    outFile = open('ConsistanceResults/results_'+subset+'_'+args.tumor+'.csv','w')
    outFile.write('\n'.join(writeMe))
    outFile.close()

