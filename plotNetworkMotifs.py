##########################################################
## Consistilator:  plotNetworkMotifs.py                 ##
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
import numpy as np
import pandas as pd
import networkx as nx
#from networkx.algorithms import community
from networkx.algorithms import clique
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from lifelines import KaplanMeierFitter
from lifelines.statistics import pairwise_logrank_test
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
import boolean2
from boolean2 import util, state, network
import palettable as pal
import pyBinarize as pyBin
import os
#import copy
#from multiprocessing import Pool, cpu_count, Manager
from matplotlib.patches import FancyArrowPatch, Circle
#import pdb
import argparse

plt.switch_backend('agg')

parser = argparse.ArgumentParser(description='use for testing argparse.py')
parser.add_argument('--tumors', help='List of tumor types', type = str)
args = parser.parse_args()

def simulation(model, trans):
    "One simulation step will update the transition graph"

    # generates all states, set limit to a value to keep only the first that many states
    # when limit is a number it will take the first that many initial states
    initializer = state.all_initial_states( model.nodes, limit=None )

    # the data is the inital data, the func is the initializer
    for data, initfunc in initializer:
        model.initialize(missing=initfunc)
        model.iterate(100)
        trans.add( model.states, times=range(100) )
    return trans

## Plot network information
##  _____________________________________
##  |        |        |        |        |
##  | NetMot |  R^2   |  And   |   Or   |
##  |        |        | Attr.  |  Attr. |
##  _____________________________________
##  |        |        | Surv.  | Surv.  |    [ TCGA LGG_GBM
##  | Dist   | Surv.  |  And   |   Or   | X3 [ REMBRANDT
##  | States | States | Attr.  | Attr.  |    [ French, et al.
##  _____________________________________
##  |        |        | Surv.  | Surv.  |    [ TCGA LGG_GBM
##  | Dist   | Surv.  |  And   |   Or   | X3 [ REMBRANDT
##  | States | States | Attr.  | Attr.  |    [ French, et al.
##  _____________________________________
#def plotNetworks(geneSets, ids, consistencies, networks, rsq, exp, binExp, pheno, symbol2entrez):
#def plotNetworks(geneSets, ids, consistencies, networks, rules_and, rules_or, rsq, exp, binExp, pheno, exp_lgg, binExp_lgg, pheno_lgg, symbol2entrez):
def plotNetworks(geneSets, ids, consistencies, networks, rules_and, rules_or, rsq, binExp_lgg, pheno_lgg, symbol2entrez):#, samples, binExp_sc, pheno_lgg, binExp_lgg
    nodeColors = ['k','r','g']
    pp = PdfPages('PipeOut/'+args.tumors+'NetMotifs_3node_seq.pdf')
    outFile = open('PipeOut/pairwise_comparisons.csv','w')
    states = ['000', '100','010','001','110','011','101','111']
    for set1 in range(len(geneSets)):
        numRows=1+len(binExp_lgg)
        fig = plt.figure(figsize=(20,(2+numRows*2.5)))
        plt.rcParams['legend.fontsize'] = 8
        grid = plt.GridSpec(numRows,5, wspace=0.25, hspace=0.25, left=0.075, right=0.95, bottom=0.05, top=0.95)
        
        # Make boolean networks first to ensure order is conserved across all aspects
        # And
        model_and = boolean2.Model( text=rules_and[set1], mode='sync')
        #trans_and = network.TransGraph( logfile='threenodes.log', verbose=True )
        trans_and = network.TransGraph( verbose=True )
        simulation( model_and, trans_and)
        att_and_graph = nx.DiGraph()
        for state_edge in trans_and.graph.adjacency():
            att_and_graph.add_edge(state_edge[0],state_edge[1].keys()[0])
            print state_edge[0],state_edge[1].keys()[0]
        attractors_and = sorted(nx.connected_components(att_and_graph.to_undirected()), key=len, reverse=True)
        
        # Or
        model_or = boolean2.Model( text=rules_or[set1], mode='sync')
        #trans_or = network.TransGraph( logfile='threenodes.log', verbose=True )
        trans_or = network.TransGraph( verbose=True )
        simulation (model_or, trans_or)
        att_or_graph = nx.DiGraph()
        for state_edge in trans_or.graph.adjacency():
            att_or_graph.add_edge(state_edge[0],state_edge[1].keys()[0])
            print state_edge[0],state_edge[1].keys()[0]
        attractors_or = sorted(nx.connected_components(att_or_graph.to_undirected()), key=len, reverse=True)
        
        if not model_and.states[0].keys()==model_or.states[0].keys():
            print geneSets[set1], model_and.states[0].keys(), model_or.states[0].keys()
            break
        nodes = dict(zip(model_and.states[0].keys(),nodeColors))

        # Network plot [0,0]
        plt.subplot(grid[0,0])
        G = networks[set1]
        node_colors = [nodes[i] for i in list(G.nodes)]
        edges = G.edges()
        edge_colors = [G[u][v]['color'] for u,v in edges]
        nx.draw(G, pos=nx.circular_layout(G),label_pos=3,with_labels=False,node_size=500,node_color=node_colors,edge_color=edge_colors,width=3,arrowsize=25,font_color='w')
        #createing label offset on network figure
        label_ratio = 0.27
        pos_labels = {} 
        #For each node in the Graph
        pos = nx.circular_layout(G)
        p=0
        for aNode in G.nodes():
            #Get the node's position from the layout
            x,y = pos[aNode]
            #Set Offset
            if p==0:
                pos_labels[aNode] = (x-label_ratio, y)
            else:
                pos_labels[aNode] = (x+label_ratio, y)
            p+=1
        nx.draw_networkx_labels(G,pos=pos_labels,fontsize=3)
        plt.title('Network Motif',fontdict={'fontsize':8})

         #Plot R squared values [0,1]
        ax = plt.subplot(grid[0,1])
        index = np.arange(len(model_and.states[0].keys()))
        print index
        plt.xticks(index, model_and.states[0].keys())
        g1, g2, g3 = plt.bar(index, [consistencies[set1]['all'][j] for j in model_and.states[0].keys()])
        print consistencies[set1]['all']
        plt.ylim((0,1))
        plt.ylabel('$R^2$')
        #plt.xlabel('Gene')
        ax.axhline(rsq,color='k',linestyle='--',alpha=0.5)
        g1.set_facecolor(nodes[list(model_and.states[0].keys())[0]])
        g2.set_facecolor(nodes[list(model_and.states[0].keys())[1]])
        g3.set_facecolor(nodes[list(model_and.states[0].keys())[2]])
        plt.title('$R^2$ for Inputs',fontdict={'fontsize':8})
        plt.tight_layout()

        # And attractors [0,2]
        print attractors_and
        plt.subplot(grid[0,2])
        plt.xticks([])
        plt.yticks([])
        positions={'000':[-1,-1],'001':[-1,0],'010':[-1,1],'100':[0,-1],'110':[0,0],'011':[0,1],'101':[1,-1],'111':[1,0]}
        pos = nx.spring_layout(att_and_graph,k=0.8,pos=positions,iterations=10)
        #pos = nx.nx_agraph.graphviz_layout(att_and_graph, prog='dot') ### TODO ### See if we can get this to work.
        selfies = list(att_and_graph.nodes_with_selfloops())
        nx.draw_networkx_nodes(att_and_graph,pos,nodelist=[i for i in att_and_graph.nodes if not i in selfies],node_color='w',node_size=100)
        nx.draw_networkx_nodes(att_and_graph,pos,nodelist=selfies,node_color='#fa9fb5',node_size=100)
        nx.draw_networkx_edges(att_and_graph,pos,fontsize=3, alpha=0.5)
        nx.draw_networkx_labels(att_and_graph,pos,fontsize=3,font_color='k')
        plt.title('AND Attractors',fontdict={'fontsize':8})
        
        # Or attractors [0,3]
        print attractors_or
        plt.subplot(grid[0,3])
        plt.xticks([])
        plt.yticks([])
        positions={'000':[-1,-1],'001':[-1,0],'010':[-1,1],'100':[0,-1],'110':[0,0],'011':[0,1],'101':[1,-1],'111':[1,0]}
        pos = nx.spring_layout(att_or_graph,k=0.8,pos=positions,iterations=10)
        selfies = list(att_or_graph.nodes_with_selfloops())
        nx.draw_networkx_nodes(att_or_graph,pos,nodelist=[i for i in att_or_graph.nodes if not i in selfies],node_color='w',node_size=100)
        nx.draw_networkx_nodes(att_or_graph,pos,nodelist=selfies,node_color='#fa9fb5',node_size=100)
        nx.draw_networkx_edges(att_or_graph,pos,fontsize=3, alpha=0.5)
        nx.draw_networkx_labels(att_or_graph,pos,fontsize=3,font_color='k')
        plt.title('OR Attractors',fontdict={'fontsize':8})
        
        # Make binary data for each
        states = ['000', '100','010','001','110','011','101','111']
        states_colors = dict(zip(states,['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray']))
        plt.tight_layout()

        ### LGG+GBM samples
        cohorts = ['TCGA']#,'GSE',REMBRANDT','French']
        for addCohorts in runFile:
            cohorts.append(addCohorts)
        header1 = '\nStudy,Gene1,Gene2,Gene3,'
        tmp1 = []
        for i in range(0,7):
            for j in range(i+1,8):
                tmp1.append(states[i]+'_vs_'+states[j])
        for i in range(0,len(attractors_and)-1):
            for j in range(i+1,len(attractors_and)):
                tmp1.append('_'.join([str(k) for k in attractors_and[i]])+'_vs_'+'_'.join([str(k) for k in attractors_and[j]])+'_and')
        for i in range(0,len(attractors_or)-1):
            for j in range(i+1,len(attractors_or)):
                tmp1.append('_'.join([str(k) for k in attractors_or[i]])+'_vs_'+'_'.join([str(k) for k in attractors_or[j]])+'_or')

        header1 = header1+','.join(tmp1)
#        print header1
        outFile.write(header1)
        for cohort in range(len(binExp_lgg)):#range(len(cohorts)):::range(len(runFile)):::len(binExp_lgg)
            print cohort, cohorts[cohort]
            # Bin for LGG+GBM
            if len([i for i in binExp_lgg[cohorts[cohort]].index if i in model_and.states[0].keys()])==3:
                tmp = binExp_lgg[cohorts[cohort]].loc[model_and.states[0].keys()].transpose()
                tmp1 = dict(zip(tmp.index,[''.join([str(int(k)) for k in tmp.loc[j].values]) for j in tmp.index]))
                tmp2 = []
                for j in pheno_lgg[cohorts[cohort].split('_')[0]].index:
                    if not j in tmp1:
                        tmp2.append(np.nan)
                    else:
                        tmp2.append(tmp1[j])
                pheno2 = pheno_lgg[cohorts[cohort].split('_')[0]].assign(binSet=tmp2)

                # Distribution of LGG+GBM states [1,0][0] as heatmap
                grades = list(pd.Series(list(set(list(pheno_lgg[cohorts[cohort].split('_')[0]]['subtype'])))).dropna())#'TCGA'->cohorts[cohort]
                #grades = ['All','IV','III','II']#replace with subtype, set of vals, exclude NAs
                #ax = plt.subplot(inner[0])
                groups = pheno2['binSet']
                dist1 = groups.value_counts()
                all1 = []
                total = float(sum(list(dist1)))
                for i in states:
                    if i in dist1:
                        all1.append(float(dist1[i])/total)
                    else:
                        all1.append(0)
                #loop through all subtypes
                gtype = {}
                for subtype in grades:
                    subset1 = pheno2.loc[pheno2['subtype']==subtype].index
                    groups = pheno2['binSet'].loc[subset1]
                    dist1 = groups.value_counts()
                    gtype[str(subtype)] = []
                    total = float(sum(list(dist1)))
                    for i in states:
                        if i in dist1:
                            gtype[str(subtype)].append(float(dist1[i])/total)
                        else:
                            gtype[str(subtype)].append(0)
                    arrayList = []
                    arrayList.append(all1)
                    for a_list in gtype:
                        arrayList.append(gtype[a_list])
                
                # Heatmap
                props1 = np.array(arrayList)
                ax = plt.subplot(grid[cohort+1,0])
                im = ax.imshow(props1, cmap=plt.get_cmap('Reds'), vmin=0, vmax=0.5, aspect='auto')
                ax.set_xticks(np.arange(len(states)))
                ax.set_yticks(np.arange(len(grades)))
                ax.set_xticklabels(states)
                ax.set_yticklabels(grades)
                plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
                for i in range(len(grades)):
                    for j in range(len(states)):
                        if props1[i,j]>0.25:
                            col1 = 'w'
                        else:
                            col1 = 'k'
                        text = ax.text(j, i, '%s' % float('%.2g' % props1[i, j]), ha="center", va="center", color=col1, fontdict={'fontsize':8})
                ax.set_title(cohorts[cohort],fontdict={'fontsize':8})

                # LGG+GBM Survival [1,1]
                ax = plt.subplot(grid[cohort+1,1])
                kmf = KaplanMeierFitter()
                pheno3 = pheno2[['OS.time','OS','binSet']].dropna(axis='rows')
                T = pheno3['OS.time']
                E = pheno3['OS']==1
                groups = pheno3['binSet']
                print 'length of groups'
                print len(groups)
                pw1 = pairwise_logrank_test(event_durations=T, event_observed=E, groups = groups)#here, divide by zero
                #print pw1.keys()
                # Add to file
                writeMe = [cohorts[cohort]] + model_and.states[0].keys()
                for state1 in range(0,7):
                    if states[state1] in pw1.keys():
                        #print pw1[states[state1]].keys()
                        for state2 in range(state1+1,8):
                            if states[state2] in pw1[states[state1]].keys():
                                #print states[state1], states[state2]
                                writeMe.append(str(pw1[states[state1]][states[state2]].p_value))
                            else:
                                writeMe.append('NA')
                    else:
                        for state2 in range(state1+1,8):
                            writeMe.append('NA')
                #print tmp
                dist1 = groups.value_counts()
                for k in states:
                    if k in dist1:
                        ix = (groups==k)
                        if len(T[ix]) > 0:
                            kmf.fit(T[ix], E[ix], label=k)
                            kmf.plot(ax=ax, ci_show=False, color=states_colors[k])
                plt.title('Survival Bin. States ('+cohorts[cohort]+')',fontdict={'fontsize':8})

                # Plot both attractors
                inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=grid[cohort+1,2], wspace=0.3, hspace=0.7)
                
                ## Attractor distributions
                # And
                attrSets = [', '.join(list(i)) for i in attractors_and]
                attrDict = dict(zip(attrSets,attractors_and))
                tmp1 = []
                for i in pheno2['binSet']:
                    for j in attrDict:
                        if i in attrDict[j]:
                            tmp1.append(j)
                attrAnd = [float(i)/float(len(tmp1)) for i in pd.Series(tmp1).value_counts()]
                for i in attrSets:# removes filtered item from corrosponding list, attrSets.
                    if i not in tmp1:
                        attrSets.remove(i)
                # Heatmap
                ax = plt.subplot(inner[0])
                props1 = np.array([attrAnd])
                print props1
                im = ax.imshow(props1, cmap=plt.get_cmap('Reds'), vmin=0, vmax=1, aspect='auto')
                ax.set_xticks(np.arange(len(attrSets)))
                ax.set_yticks(np.arange(1))
                ax.set_xticklabels(attrSets)
                ax.set_yticklabels(['And'])
                plt.setp(ax.get_xticklabels(), rotation=15, ha="center", rotation_mode="anchor")
                for i in range(1):
                    for j in range(len(attrSets)):#quick fix to run code, added -1, due to out of index
                        if props1[i,j]>0.7:
                            col1 = 'w'
                        else:
                            col1 = 'k'
                        text = ax.text(j, i, '%s' % float('%.2g' % props1[i, j]), ha="center", va="center", color=col1, fontdict={'fontsize':8})
                ax.set_title(cohorts[cohort],fontdict={'fontsize':8})

                # Or
                attrSets = [', '.join(list(i)) for i in attractors_or]
                attrDict = dict(zip(attrSets,attractors_or))
                tmp1 = []
                for i in pheno2['binSet']:
                    for j in attrDict:
                        if i in attrDict[j]:
                            tmp1.append(j)
                attrOr = [float(i)/float(len(tmp1)) for i in pd.Series(tmp1).value_counts()]
                for i in attrSets:# removes filtered item from corrosponding list, attrSets.
                    if i not in tmp1:
                        attrSets.remove(i)
                # Heatmap
                ax = plt.subplot(inner[1])
                props1 = np.array([attrOr])
                print props1
                im = ax.imshow(props1, cmap=plt.get_cmap('Reds'), vmin=0, vmax=1, aspect='auto')
                ax.set_xticks(np.arange(len(attrSets)))
                ax.set_yticks(np.arange(1))
                ax.set_xticklabels(attrSets)
                ax.set_yticklabels(['Or'])
                plt.setp(ax.get_xticklabels(), rotation=15, ha="center", rotation_mode="anchor")
                for i in range(1):
                    for j in range(len(attrSets)):#quick fix to run code, added -1, due to out of index
                        if props1[i,j]>0.7:
                            col1 = 'w'
                        else:
                            col1 = 'k'
                        text = ax.text(j, i, '%s' % float('%.2g' % props1[i, j]), ha="center", va="center", color=col1, fontdict={'fontsize':8})
                #ax.set_title(cohorts[cohort],fontdict={'fontsize':8})

                # LGG+GBM Survival grouped by AND attractor states [1,1]
                ax = plt.subplot(grid[cohort+1,3])
                kmf = KaplanMeierFitter()
                pheno3 = pheno2[['OS.time','OS','binSet']].dropna(axis='rows')
                T = pheno3['OS.time']
                E = pheno3['OS']==1
                groups = pheno3['binSet']
                dist1 = groups.value_counts()
                attrSets = [', '.join(list(i)) for i in attractors_and]
                attrDict = dict(zip(attrSets,attractors_and))
                #print attrSets
                #print attrDict
                attrColors = dict(zip(attrSets,pal.tableau.Tableau_10.hex_colors[0:len(attractors_and)]))
                groupsAttr = []
                for stater in groups:
                    for k in attrDict:
                        if stater in attrDict[k]:
                            groupsAttr.append(k)
                #print groupsAttr
                print 'length of groupsAttr'
                print len(groupsAttr)
#                for echGroups in range(len(groupsAttr)):
#                    print echGroups
#                print 'END'
                pw1 = pairwise_logrank_test(event_durations=T, event_observed=E, groups = groups)#here broken
                # Add to file
                #writeMe = [cohorts[cohort]] + model_and.states[0].keys()
                for state1 in range(0,(len(attrSets)-1)):
                    if attrSets[state1] in pw1.keys():
                        for state2 in range(state1+1,len(attrSets)):
                            if attrSets[state2] in pw1[attrSets[state1]].keys():
                                writeMe.append(str(pw1[attrSets[state1]][attrSets[state2]].p_value))
                            else:
                                writeMe.append('NA')
                    else:
                        for state2 in range(state1+1,len(attrSets)):
                            writeMe.append('NA')
                for k in attrDict:
                    ix = []
                    for i in groups:
                        if not i in attrDict[k]:
                            ix.append(False)
                        else:
                            ix.append(True)
                    #print ix
                    #print k, ix, T[ix], E[ix]
                    if len(T[ix])>0:
                        kmf.fit(T[ix], E[ix], label=k)
                        kmf.plot(ax=ax, ci_show=False, color=attrColors[k])
                plt.title('Survival AND Attractor States ('+cohorts[cohort]+')',fontdict={'fontsize':8})
                #print dir(ax)
                #print ax._get_lines().get_next_color()
                
                # LGG+GBM Survival grouped by OR attractor states [1,1]
                ax = plt.subplot(grid[cohort+1,4])
                kmf = KaplanMeierFitter()
                pheno3 = pheno2[['OS.time','OS','binSet']].dropna(axis='rows')
                T = pheno3['OS.time']
                E = pheno3['OS']==1
                groups = pheno3['binSet']
                dist1 = groups.value_counts()
                attrSets = [', '.join(list(i)) for i in attractors_or]
                attrDict = dict(zip(attrSets,attractors_or))
                attrColors = dict(zip(attrSets,pal.tableau.Tableau_10.hex_colors[0:len(attractors_or)]))
                groupsAttr = []
                for stater in groups:
                    for k in attrDict:
                        if stater in attrDict[k]:
                            groupsAttr.append(k)
                pw1 = pairwise_logrank_test(event_durations=T, event_observed=E, groups = groupsAttr)
                # Add to file
                #writeMe = [cohorts[cohort]] + model_and.states[0].keys()
                for state1 in range(0,(len(attrSets)-1)):
                    if attrSets[state1] in pw1.keys():
                        for state2 in range(state1+1,len(attrSets)):
                            if attrSets[state2] in pw1[attrSets[state1]].keys():
                                writeMe.append(str(pw1[attrSets[state1]][attrSets[state2]].p_value))
                            else:
                                writeMe.append('NA')
                    else:
                        for state2 in range(state1+1,len(attrSets)):
                            writeMe.append('NA')
                for k in attrDict:
                    ix = []
                    for i in groups:
                        if not i in attrDict[k]:
                            ix.append(False)
                        else:
                            ix.append(True)
                    #print ix
                    #print k, ix, T[ix], E[ix]
                    if len(T[ix])>0:
                        kmf.fit(T[ix], E[ix], label=k)
                        kmf.plot(ax=ax, ci_show=False, color=attrColors[k])
                plt.title('Survival OR Attractor States ('+cohorts[cohort]+')',fontdict={'fontsize':8})
                
                # Write out survival results
                outFile.write('\n'+','.join(writeMe))
                

        ## Plot single cell data
#        scData = ['GSE57872','GSE84465','GSE89567','GSE102130','GSE70630']
#        
#        scData = ['IV','III','II','H3K27M']
#        for dataset in range(len(scData)):
#            sources = sorted(list(set(samples[scData[dataset]])))
#            print len(sources), sources
#            
#            # Concatenate states
#            tmp = binExp_sc[scData[dataset]].loc[model_and.states[0].keys()].transpose().dropna()
#            inc_samps = [i for i in range(len(binExp_sc[scData[dataset]].columns)) if binExp_sc[scData[dataset]].columns[i] in tmp.index]
#            #print scData[dataset], inc_samps
#            tmp_samples = [samples[scData[dataset]][i] for i in inc_samps]
#            if tmp.shape[0]>0:
#                tmp1 = tmp.assign(binSet=[''.join([str(int(k)) for k in tmp.loc[j].values]) for j in tmp.index])
#             
#                # Distribution of LGG+GBM states [1,0][0]
#                allLines = []
#                allLines_and = []
#                allLines_or = []
#                for s1 in range(len(sources)):
#                    cur = [i for i in range(len(tmp_samples)) if tmp_samples[i]==sources[s1]]
#                    if len(cur)>0:
#                        #print cur, tmp1['binSet'].shape
#                        groups = tmp1['binSet'].iloc[cur]
#                        dist1 = groups.value_counts()
#                        dist2 = []
#                        total = float(sum(list(dist1)))
#                        for i in states:
#                            if i in dist1:
#                                dist2.append(float(dist1[i])/total)
#                            else:
#                                dist2.append(0)
#                        allLines.append(dist2)
#                    
#                        # And
#                        attrSets_and = [', '.join(list(i)) for i in attractors_and]
#                        attrDict = dict(zip(attrSets_and,attractors_and))
#                        tmp2 = []
#                        for i in tmp1['binSet'].iloc[cur]:
#                            for j in attrDict:
#                                if i in attrDict[j]:
#                                    tmp2.append(j)
#                        attrAnd = []
#                        dist1 = pd.Series(tmp2).value_counts()
#                        total = float(sum(list(dist1)))
#                        for i in attrSets_and:
#                            if i in dist1:
#                                attrAnd.append(float(dist1[i])/total)
#                            else:
#                                attrAnd.append(0)
#                        allLines_and.append(attrAnd)
#
#                        # Or
#                        attrSets_or = [', '.join(list(i)) for i in attractors_or]
#                        attrDict = dict(zip(attrSets_or,attractors_or))
#                        tmp2 = []
#                        for i in tmp1['binSet'].iloc[cur]:
#                            for j in attrDict:
#                                if i in attrDict[j]:
#                                    tmp2.append(j)
#                        attrOr = []
#                        dist1 = pd.Series(tmp2).value_counts()
#                        total = float(sum(list(dist1)))
#                        for i in attrSets_or:
#                            if i in dist1:
#                                attrOr.append(float(dist1[i])/total)
#                            else:
#                                attrOr.append(0)
#                        allLines_or.append(attrOr)
#                    else:
#                        allLines.append([np.NaN,np.NaN,np.NaN,np.NaN,np.NaN,np.NaN,np.NaN,np.NaN])
#                        allLines_and.append([np.NaN for i in range(len(attractors_and))])
#                        allLines_or.append([np.NaN for i in range(len(attractors_or))])
#                
#                # Heatmap
#                props1 = np.array(allLines)
#                print(props1)
#                ax = plt.subplot(grid[4,dataset])
#                masked_array = np.ma.array(allLines, mask=np.isnan(allLines))
#                print masked_array
#                cmap = plt.get_cmap('Reds')
#                cmap.set_bad('white',1.)
#                im = ax.imshow(masked_array, interpolation='nearest', cmap=cmap, vmin=0, vmax=0.5, aspect='auto')
#                #im = ax.imshow(props1, cmap=plt.get_cmap('Reds'))
#                ax.set_xticks(np.arange(len(states)))
#                ax.set_yticks(np.arange(len(sources)))
#                ax.set_xticklabels(states)
#                ax.set_yticklabels(sources)
#                plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
#                for i in range(len(sources)):
#                    for j in range(len(states)):
#                        if props1[i,j]>0.25:
#                            col1 = 'w'
#                        else:
#                            col1 = 'k'
#                        text = ax.text(j, i, '%s' % float('%.1g' % props1[i, j]), ha="center", va="center", color=col1, fontdict={'fontsize':8})
#                ax.set_title(scData[dataset],fontdict={'fontsize':8})
#
#                ## Attractor distributions
#                # Heatmap
#                print allLines
#                print allLines_and
#                props1 = np.array(allLines_and)
#                #print props1
#                ax = plt.subplot(grid[5,dataset])
#                #print allLines_and
#                masked_array = np.ma.array (allLines_and, mask=np.isnan(allLines_and))
#                cmap = plt.get_cmap('Reds')
#                cmap.set_bad('white',1.)
#                im = ax.imshow(masked_array, interpolation='nearest', cmap=cmap, vmin=0, vmax=0.5, aspect='auto')
#                #im = ax.imshow(props1, cmap=plt.get_cmap('Reds'), vmin=0, vmax=1, aspect='auto')
#                ax.set_xticks(np.arange(len(attrSets_and)))
#                ax.set_yticks(np.arange(len(sources)))
#                ax.set_xticklabels(attrSets_and)
#                ax.set_yticklabels(sources)
#                plt.setp(ax.get_xticklabels(), rotation=15, ha="center", rotation_mode="anchor")
#                for i in range(len(sources)):
#                    for j in range(len(attrSets_and)):
#                        if props1[i,j]>0.45:
#                            col1 = 'w'
#                        else:
#                            col1 = 'k'
#                        text = ax.text(j, i, '%s' % float('%.2g' % props1[i, j]), ha="center", va="center", color=col1, fontdict={'fontsize':8})
#                ax.set_title(scData[dataset]+' AND',fontdict={'fontsize':8})
#
#                # Heatmap
#                props1 = np.array(allLines_or)
#                #print props1
#                ax = plt.subplot(grid[6,dataset])
#                masked_array = np.ma.array (allLines_or, mask=np.isnan(allLines_or))
#                cmap = plt.get_cmap('Reds')
#                cmap.set_bad('white',1.)
#                im = ax.imshow(masked_array, interpolation='nearest', cmap=cmap, vmin=0, vmax=0.5, aspect='auto')
#                #im = ax.imshow(props1, cmap=plt.get_cmap('Reds'), vmin=0, vmax=1, aspect='auto')
#                ax.set_xticks(np.arange(len(attrSets_or)))
#                ax.set_yticks(np.arange(len(sources)))
#                ax.set_xticklabels(attrSets_or)
#                ax.set_yticklabels(sources)
#                plt.setp(ax.get_xticklabels(), rotation=15, ha="center", rotation_mode="anchor")
#                for i in range(len(sources)):
#                    for j in range(len(attrSets_or)):
#                        #print props1.shape, i, j
#                        if props1[i,j]>0.45:
#                            col1 = 'w'
#                        else:
#                            col1 = 'k'
#                        text = ax.text(j, i, '%s' % float('%.2g' % props1[i, j]), ha="center", va="center", color=col1, fontdict={'fontsize':8})
#                ax.set_title(scData[dataset]+' OR',fontdict={'fontsize':8})

        if set1==6:
            plt.show()
            break
        pp.savefig(fig)
        #plt.show()
        #break
    pp.close()
    outFile.close()

# Read in Biotapesty file
inEdges = {}
inFile = open('PanOutput/biotapestry_CHIR_curve'+args.tumors+'.csv','r')
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
    for i in node2.split(';'):
        if not i in inEdges:
            inEdges[i] = {}
        if not node1 in inEdges[i]:
            inEdges[i][node1] = splitUp[6].strip('"')
inFile.close()

# Load up genesets and netMatrices
data = {}
consistent = {}
gene2probe = {}
subsets = ['all']
rsq = 0.5
for subset in subsets:
    inFile = open('ConsistanceResults/results_'+subset+'_'+args.tumors+'.csv','r')#'ConsistanceResults/results_'+subset+'_'+args.tumors+'.csv'
    inFile.readline() # Get rid of header
    while 1:
        line = inFile.readline()
        if not line:
            break
        splitUp = line.strip().split(',')
        if not (splitUp[1] in ['6', '12','36']):
            if not splitUp[1] in data:
                data[splitUp[1]] = {}
                consistent[splitUp[1]] = {}
            if not splitUp[3] in data[splitUp[1]]:
                data[splitUp[1]][splitUp[3]] = {}
                consistent[splitUp[1]][splitUp[3]] = {}
            data[splitUp[1]][splitUp[3]][subset] = splitUp
            if not (splitUp[6]=='Inconsistent' or splitUp[9]=='Inconsistent' or splitUp[12]=='Inconsistent') and (splitUp[6]=='NA' or float(splitUp[6])>=rsq) and (splitUp[9]=='NA' or float(splitUp[9])>=rsq) and (splitUp[12]=='NA' or float(splitUp[12])>=rsq):
                consistent[splitUp[1]][splitUp[3]][subset] = 'Yes'
            else:
                consistent[splitUp[1]][splitUp[3]][subset] = 'No'
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

# Read in binarized data args.tumors
binExp_lgg = {}
binExp_lgg['TCGA'] = pd.read_csv('exprs_all/'+args.tumors+'_RNAseq.csv', header=0, index_col=0)#.transpose()
binExp_lgg['TCGA'].columns = [i.replace('.','-') for i in binExp_lgg['TCGA'].columns.values]
binExp_lgg['TCGA']=binExp_lgg['TCGA'].transpose().dropna().transpose()
binExp_lgg['TCGA'] = pyBin.binarize_kMeans_matrix(binExp_lgg['TCGA'])
binIndex = binExp_lgg['TCGA'].index.values
newIndex = []
count = 0
for geneid in binIndex:
    if entrez2symbol.has_key(str(geneid)):
        newIndex.append(entrez2symbol[str(geneid)])
    else:
        newIndex.append(geneid)
    count = count+1
#binExp_lgg['TCGA'] = pd.read_csv('ACC_TCGA_Binarize.csv', header=0, index_col=0)
binExp_lgg['TCGA'].index=newIndex
#binExp_lgg['TCGA'].to_csv('ACC_TCGA_Binarize.csv', index=False)

#Examine expr files in dir
for path, subdirs, files in os.walk('/home/ssstrike/Documents/TF_Scope/Reconstructed_exprs'):#C:\Users\Shawn\Documents\Research\PanCancer\Reconstructed_exprs:::/home/ssstrike/Documents/TF_Scope/Reconstructed_exprs
    print files
runFile=[]
for fileName in files:
    splitName=fileName.split('_')
    if splitName[0] == args.tumors:#args.tumors
        runFile.append(splitName[1]+'_'+splitName[2])

#loop for modular number of key inputs
for numExpFile in range(len(runFile)):
    runName=runFile[numExpFile].split('_')[0]+'_'+runFile[numExpFile].split('_')[1]
    binExp_lgg[runName] = pd.read_csv('Reconstructed_exprs/'+args.tumors+'_'+runFile[numExpFile].split('_')[0]+'_'+runFile[numExpFile].split('_')[1]+'_exprs_entrez.csv', header=0, index_col=0)#.transpose()
    binExp_lgg[runName].columns = [i.replace('.','-') for i in binExp_lgg[runName].columns.values]
    binExp_lgg[runName]=binExp_lgg[runName].transpose().dropna().transpose()
    binExp_lgg[runName] = pyBin.binarize_kMeans_matrix(binExp_lgg[runName])
    binIndex = binExp_lgg[runName].index.values
    newIndex = []
    count = 0
    for geneid in binIndex:
        if entrez2symbol.has_key(str(geneid)):
            newIndex.append(entrez2symbol[str(geneid)])
        else:
            newIndex.append(geneid)
        count = count+1
#    binExp_lgg[runName] = pd.read_csv(args.tumors+'_'+runFile[numExpFile].split('_')[0]+'_'+runFile[numExpFile].split('_')[1]+'_Binarize.csv', header=0, index_col=0)
    binExp_lgg[runName].index=newIndex
#    binExp_lgg[runName].to_csv(args.tumors+'_'+runFile[numExpFile].split('_')[0]+'_'+runFile[numExpFile].split('_')[1]+'_Binarize.csv', index=False)


#binExp_lgg['REMBRANDT'] = pd.read_csv('../gbm_lgg/binGeneExp_REMBRANDT_VMedian.csv', header=0, index_col=0)#.transpose()
#binExp_lgg['REMBRANDT'] = binExp_lgg['REMBRANDT'].dropna()
#binExp_lgg['French'] = pd.read_csv('../gbm_lgg/binGeneExp_French_VMedian.csv', header=0, index_col=0)#.transpose()
#binExp_lgg['French'] = binExp_lgg['French'].dropna()att


# Phenotypes
pheno_lgg = {}
pheno_lgg['TCGA'] = pd.read_csv('phenotypes_panCancer.csv', header=0, index_col=0)
pheno_lgg['TCGA'] = pheno_lgg['TCGA'].loc[pheno_lgg['TCGA']['tumor']==args.tumors]#args.tumors
pheno_lgg['TCGA'] = pheno_lgg[pheno_lgg.keys()[0]][["tumor","OS.time","OS","subtype"]]

pheno_lgg['GSE'] = pd.read_csv('GSEPheno_panCancer.csv', header=0, index_col=0)
pheno_lgg['GSE'] = pheno_lgg['GSE'].loc[pheno_lgg['GSE']['tumor']==args.tumors]
pheno_lgg['GSE'] = pheno_lgg[pheno_lgg.keys()[1]][["tumor","OS.time","OS","subtype"]]
    
#pheno_lgg['REMBRANDT'] = pd.read_csv('../gbm_lgg/phenotypes_REMBRANDT.csv', header=0, index_col=0)
#pheno_lgg['French'] = pd.read_csv('../gbm_lgg/phenotypes_French.csv', header=0, index_col=0)

# Read in single cell RNA-seq data
# Read in single cell RNA-seq data

#binExp_sc = {}
#samples = {}
#tmp = pd.read_csv('../singleCell/binGeneExp_scGBM_Vgt0_GSE57872.csv', header=0, index_col=0)#.transpose()
#binExp_sc['IV'] = tmp
#samples['IV'] = [i.split('_')[0] for i in tmp.columns]
#tmp = pd.read_csv('../singleCell/binGeneExp_scGBM_Vgt0_GSE84465.csv', header=0, index_col=0)#.transpose()
#binExp_sc['IV'] = pd.concat([binExp_sc['IV'],tmp], axis=1)
#samples['IV'] = samples['IV']+['BT_S2']*1169+['BT_S1']*489+['BT_S4']*1542+['BT_S6']*389
#tmp = pd.read_csv('../singleCell/binGeneExp_scGBM_Vgt0_GSE89567.csv', header=0, index_col=0)#.transpose()
#tmp2 = [j.split('.')[0] for j in [i.split('_')[0] for i in tmp.columns]]
#for i in range(len(tmp2)):
#    if tmp2[i]=='X57':
#        tmp2[i] = 'MGH57'
#    if tmp2[i]=='mgh103':
#        tmp2[i] = 'MGH103'
#gbmIV = [i for i in range(len(tmp2)) if tmp2[i] in ['MGH45','MGH57']]
#binExp_sc['IV'] = pd.concat([binExp_sc['IV'],tmp[tmp.columns[gbmIV]]], axis=1)
#samples['IV'] = samples['IV']+[tmp2[i] for i in gbmIV]
#gliomaIII = [i for i in range(len(tmp2)) if tmp2[i] in ['MGH42','MGH43','MGH44','MGH56','MGH61','MGH64','MGH103']]
#binExp_sc['III'] = tmp[tmp.columns[gliomaIII]]
#samples['III'] = [tmp2[i] for i in gliomaIII]
#gliomaII = [i for i in range(len(tmp2)) if tmp2[i] in ['MGH107neg','MGH107pos']]
#binExp_sc['II'] = tmp[tmp.columns[gliomaII]]
#samples['II'] = [tmp2[i] for i in gliomaII]
#binExp_sc['H3K27M'] = pd.read_csv('../singleCell/binGeneExp_scGBM_K27M_Vgt0_GSE102130.csv', header=0, index_col=0)#.transpose()
#samples['H3K27M'] = [i.split('_')[0].split('.')[0] for i in binExp_sc['H3K27M'].columns]
#tmp = pd.read_csv('../singleCell/binGeneExp_scOligo_Vgt0_GSE70630.csv', header=0, index_col=0)#.transpose()
#tmp2 = [i.split('_')[0] for i in tmp.columns]
#for i in range(len(tmp2)):
#    if tmp2[i]=='X93':
#        tmp2[i] = 'MGH93'
#    if tmp2[i]=='X97':
#        tmp2[i] = 'MGH97'
#gliomaII = [i for i in range(len(tmp2)) if tmp2[i] in ['MGH36','MGH53','MGH54','MGH60','MGH93','MGH97']]
#binExp_sc['II'] = pd.concat([binExp_sc['II'],tmp[tmp.columns[gliomaII]]], axis=1)
#samples['II'] = samples['II']+[tmp2[i] for i in gliomaII]

# Gather data
geneSets = []
ids = []
netMatrices = []
consistencies = []
networks = []
rules_and = []
rules_or = []
for netMotif in data:
    for instance in data[netMotif]:
        # Only plot if significant amount of variance explained
        if len([i for i in subsets if consistent[netMotif][instance][i]=='Yes']) > 0:
            genes = data[netMotif][instance]['all'][3].split(';')
            if len(set(genes))==len(genes):
                #if len([i for i in genes if i in ['MYB','SMAD9','SPDEF']])==3:
                geneSets.append(genes)
                ids.append('id'+data[netMotif][instance]['all'][1]+' '+data[netMotif][instance]['all'][2])
                tmp = {'all':dict(zip(genes, [data[netMotif][instance]['all'][i] for i in [6,9,12]]))}
                for i in tmp:
                    for j in tmp[i]:
                        if tmp[i][j]=='Inconsistent':
                            tmp[i][j] = -0.25
                        elif not tmp[i][j]=='NA':
                            tmp[i][j] = float(tmp[i][j])
                        else:
                            tmp[i][j] = 0
                consistencies.append(tmp)

                # Make networks
                G = nx.DiGraph()
                tmp_and = ""
                tmp_or = ""
                for gene1 in genes:
                    tmp_pos = []
                    tmp_neg = []
                    if gene1 in inEdges and len([i for i in inEdges[gene1] if (i in genes and not i==gene1)])>0:
                        for gene2 in inEdges[gene1]:
                            if gene2 in genes and not gene2==gene1:
                                if inEdges[gene1][gene2]=='positive':
                                    G.add_edge(gene2,gene1,color='g')
                                    tmp_pos.append(str(gene2))
                                elif inEdges[gene1][gene2]=='negative':
                                    G.add_edge(gene2,gene1,color='r')
                                    tmp_neg.append(str(gene2))
                    else:
                        print gene1, genes
                        tmp_pos.append(gene1)
                    if not len(tmp)==0:
                        tmp_neg2_and = []
                        tmp_neg2_or = []
                        if len(tmp_neg)>0:
                            tmp_neg2_and = ['((not ('+' and '.join(tmp_neg)+')) and '+gene1+')']
                            tmp_neg2_or = ['((not ('+' or '.join(tmp_neg)+')) and '+gene1+')']
                        tmp_and += '\n'+str(gene1)+'* = '+' and '.join(tmp_pos+tmp_neg2_and)
                        tmp_or += '\n'+str(gene1)+'* = '+' or '.join(tmp_pos+tmp_neg2_or)
                rules_and.append(tmp_and)
                print tmp_and
                rules_or.append(tmp_or)
                print tmp_or
                networks.append(G)
                
# Plot them
plotNetworks(geneSets, ids, consistencies, networks, rules_and, rules_or, rsq, binExp_lgg, pheno_lgg, symbol2entrez)#, samples, binExp_sc

