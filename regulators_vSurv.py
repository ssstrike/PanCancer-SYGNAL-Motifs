import copy
import argparse
from scipy.stats import hypergeom

parser = argparse.ArgumentParser(description='use for testing argparse.py')

parser.add_argument('--correlation', help='Run at a specific correlation', type = float)
parser.add_argument('--tumors', help='List of tumor types', type = str)
args = parser.parse_args()

tumors = [item for item in args.tumors.split(',')]

#tumors = [item for item in tumors.split(',')]

#tumors = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM'] # 'GBM',
#tumors = args.tumors

runs = ['pita', 'targetscan', 'tfbs_db']

# Load up included GBM biclusters
with open('gbm/gbm_includedBiclusters.csv','r') as inFile:
    gbmInclude = [i.strip() for i in inFile.readlines()]

# Load TF conversion dictionary
tfConv = {}
with open('id_conversion/gene2entrez.csv','r') as inFile:
    while 1:
        line = inFile.readline()
        if not line:
            break
        splitUp = line.strip().split(',')
        tfConv[splitUp[0]] = splitUp[1]

# Load miRNA conversion dictionary
miRNAConv = {}
with open('gbm/mature_hsa.csv','r') as inFile:
    while 1:
        line = inFile.readline()
        if not line:
            break
        splitUp = line.strip().split(',')
        miRNAConv[splitUp[0]] = splitUp[1]

# Cutoffs
varExp = args.correlation
alphaFPC = 0
alphaMTC = varExp/63354
alphaSurv = varExp

#infiltrates = ["B.cells.naive", "B.cells.memory", "Plasma.cells", "T.cells.CD8", "T.cells.CD4.naive", "T.cells.CD4.memory.resting", "T.cells.CD4.memory.activated", "T.cells.follicular.helper", "T.cells.regulatory..Tregs.", "T.cells.gamma.delta", "NK.cells.resting", "NK.cells.activated", "Monocytes", "Macrophages.M0", "Macrophages.M1", "Macrophages.M2", "Dendritic.cells.resting", "Dendritic.cells.activated", "Mast.cells.resting", "Mast.cells.activated", "Eosinophils", "Neutrophils", "TotalLeukocyte", "OS", "OS.covAge","DSS", "DSS.covAge","PFI", "PFI.covAge"]

print 'Load in genes...'
# Load each postProcessed file
genes = {}
tumorGenes = {}
for tumor in tumors:
    tumorGenes[tumor] = []
    print tumor
    with open('postProcessed_vSurv/postProcessed_'+tumor+'_pita.csv','r') as inFile:
        header = inFile.readline().strip().split(',') # Header
        while 1:
            inLine = inFile.readline()
            if not inLine:
                break
            tmp1 = inLine.strip().split(',')
            tmp = dict(zip(header, [i.strip('"') for i in tmp1]))
            tumorGenes[tumor] += tmp['Genes'].split(' ')
    tumorGenes[tumor] = list(set(tumorGenes[tumor]))

# Load each postProcessed file
#retain = []
biclusters = {}
regulators = {}
tumorRegulators = {}
qualityBics = []
allBics = []
header = ''
btn3aWrite =[]
immmodBicsWrite =[]
for tumor in tumors:
    for run in runs:
        print tumor, run
        with open('postProcessed_vSurv/postProcessed_'+tumor+'_'+run+'.csv','r') as inFile:
            header = inFile.readline().strip().split(',') # Header
            while 1:
                inLine = inFile.readline()
                if not inLine:
                    break
                tmp1 = inLine.strip().split(',')
                tmp = dict(zip(header, [i.strip('"') for i in tmp1]))
                tmp['tumor'] = tumor
                tmp['run'] = run
                genes = tmp['Genes'].split(' ')
                # Filter for those to retain
                allBics.append(tumor+'_'+run+'_'+tmp['Bicluster'])
                if tumor=='GBM':
                    if not run+'_'+tmp['Bicluster'].strip('"') in gbmInclude:
                        continue

                if (not tmp[tumor+' Var. Exp. First PC']=='NA') and float(tmp[tumor+' Var. Exp. First PC'])>=varExp and (not tmp[tumor+' Var. Exp. First PC Perm. P-Value']=='NA') and float(tmp[tumor+' Var. Exp. First PC Perm. P-Value'])>=alphaFPC and ('OS.covAgeSex_'+tumor+'.p' in tmp) and (not tmp['OS.covAgeSex_'+tumor+'.p']=='NA') and float(tmp['OS.covAgeSex_'+tumor+'.p'])>=alphaSurv: # and ((not tmp['TotalLeukocyte_'+tumor+'.p']=='NA') and float(tmp['TotalLeukocyte_'+tumor+'.p'])<=alphaMTC):
                    # TF regulators
                    tfs = []
                    if not tmp['Up.MEME Motif1 Correlated Matches_'+tumor]=='NA':
                        tfs += [i.split(':')[0] for i in tmp['Up.MEME Motif1 Correlated Matches_'+tumor].strip().split(' ')]
                    if not tmp['Up.MEME Motif2 Correlated Matches_'+tumor]=='NA':
                        tfs += [i.split(':')[0] for i in tmp['Up.MEME Motif2 Correlated Matches_'+tumor].strip().split(' ')]
                    if not tmp['Up.WEEDER Motif1 Correlated Matches_'+tumor]=='NA':
                        tfs += [i.split(':')[0] for i in tmp['Up.WEEDER Motif1 Correlated Matches_'+tumor].strip().split(' ')]
                    if not tmp['Up.WEEDER Motif2 Correlated Matches_'+tumor]=='NA':
                        tfs += [i.split(':')[0] for i in tmp['Up.WEEDER Motif2 Correlated Matches_'+tumor].strip().split(' ')]
                    if not tmp['TFBS_DB.Correlated Matches_'+tumor]=='NA':
                        tfs += [i.split(':')[0] for i in tmp['TFBS_DB.Correlated Matches_'+tumor].strip().split(' ')]
                    if tumor=='GBM':
                        tfs = [tfConv[tf] for tf in tfs if tf in tfConv]

                    # miRNAs regulators
                    miRNAs = []
                    if not tmp['3pUTR.WEEDER Motif1 Matches']=='NA':
                        miRNAs += [i.split(':')[0] for i in tmp['3pUTR.WEEDER Motif1 Matches'].strip().split(' ')]
                    if not tmp['3pUTR.WEEDER Motif2 Matches']=='NA':
                        miRNAs += [i.split(':')[0] for i in tmp['3pUTR.WEEDER Motif2 Matches'].strip().split(' ')]
                    if not tmp['3pUTR_pita.miRNAs']=='NA' and (float(tmp['3pUTR_pita.percTargets'].strip().split(' ')[0])>=0.1 and float(tmp['3pUTR_pita.pValue'])<=0.05):
                        miRNAs += [i.split(':')[0] for i in tmp['3pUTR_pita.miRNAs'].strip().split(' ')]
                    if not tmp['3pUTR_targetScan.miRNAs']=='NA' and (float(tmp['3pUTR_targetScan.percTargets'].strip().split(' ')[0])>=0.1 and float(tmp['3pUTR_targetScan.pValue'])<=0.05):
                        miRNAs += [i.split(':')[0] for i in tmp['3pUTR_targetScan.miRNAs'].strip().split(' ')]
                    if tumor=='GBM':
                        miRNAs = [miRNAConv[miRNA] for miRNA in miRNAs if miRNA in miRNAConv]

                    if not tumor in regulators:
                        regulators[tumor] = {}
                    regulators[tumor][run+'_'+tmp['Bicluster']] = list(set(tfs))+list(set(miRNAs))
                    biclusters[tumor+'_'+run+'_'+tmp['Bicluster']] = list(set(tfs))+list(set(miRNAs))
    if tumor in regulators:
        tumorRegulators[tumor] = list(set(sum(regulators[tumor].values(),[])))
        print len(tumorRegulators[tumor])


# Intersect of all tumors: Not all tumors
intTfs = set(tumorRegulators['ACC'])
for tumor in tumors:
    if tumor in tumorRegulators:
        intTfs = intTfs.intersection(tumorRegulators[tumor])

# Intersect of all tumors: Not all tumors
cntTfs = {}
for tumor in tumors:
    if tumor in tumorRegulators:
        for tf1 in tumorRegulators[tumor]:
            if not tf1 in cntTfs:
                cntTfs[tf1] = 1
            else:
                cntTfs[tf1] += 1

intTfs = []
for tf1 in cntTfs:
    if cntTfs[tf1]>=15:
        intTfs.append(tf1)

# Number of tumors TF discovered for Leukocyte fraction
with open('RegulatorsOutput/tfTumors_vSurv.csv','w') as outFile:
    outFile.write('\n'.join([str(i)+','+str(cntTfs[i]) for i in cntTfs.keys()]))

# Biclusters from intTfs
keepers = ['tf,bicluster']
for tf in intTfs:
    for bicluster in biclusters:
        if tf in biclusters[bicluster]:
            keepers.append(tf+','+bicluster)

with open('RegulatorsOutput/output_vSurv.csv','w') as outFile:
    outFile.write('\n'.join(keepers))


