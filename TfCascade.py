import json
import pandas as pd
from scipy.stats import pearsonr
import argparse

parser = argparse.ArgumentParser(description='use for testing argparse.py')

parser.add_argument('--correlation', help='Run at a specific correlation', type = float)
parser.add_argument('--family', help='Type I Family Expansion = expand to all possible TFtargs via all motifs from family members\nType II Family Expansion = If no motif for TFreg then expand to family members', type = float)
args = parser.parse_args()

# Modifies the regulators_vSurv output files to work with TfCascade
cluster2tf = {}
inputFileName = 'RegulatorsOutput/output_vSurv.csv'
with open(inputFileName,'r') as infile:
    print 'open file: '+inputFileName
    header = infile.readline() # header
    while 1:
        line = infile.readline()
        if not line:
            break
        split = line.strip().split(',')
        cancerType = ''
        counter = 0
        while 1:
            nxtChar = split[1][counter]
            if '_' in nxtChar:
                    break
            cancerType += nxtChar
            counter += 1
        
        if not cancerType in cluster2tf:
            cluster2tf[cancerType] = []
        if not split[0] in cluster2tf[cancerType]:
            cluster2tf[cancerType].append(split[0])
            
with open ('input_tf_list.csv','w') as outfile:
    outfile.write('Tumor Type,TF List')
    for tumorType in cluster2tf.keys():
        listTf = []
        for singleTf in cluster2tf[tumorType]:
            listTf.append(singleTf)
        outfile.write('\n'+tumorType+','+','.join(listTf))
print 'Done, Input List'

# Make a Biotapestry CSV file
def biotapestry(filename, data, regions):
    writeMe = []
    writeMe.append('"# Model Commands",,,,,,,,,,')
    writeMe.append('"# Command Type","Model Name","Parent Model",,,,,,,,')
    writeMe.append('"model","root",,,,,,,,,')
    for i in data:
        writeMe.append('"model","'+i+'","root",,,,,,,,')
    writeMe.append(',,,,,,,,,,')
    writeMe.append('"# Region Commands",,,,,,,,,,')
    writeMe.append('"# Command Type","Model Name","Region Name","Region Abbreviation",,,,,,,')
    for i in data:
        writeMe.append('"region","'+i+'","A","A",,,,,,,')
        #for j in regions:
        #    writeMe.append('"region","'+i+'","Cluster-'+j+'","'+j+'",,,,,,,')
    writeMe.append(',,,,,,,,,,')
    writeMe.append('"# Standard Interactions",,,,,,,,,,')
    writeMe.append('"# Command Type","Model Name","Source Type","Source Name","Target Type","Target Name","Sign","Source Region Abbrev","Target Region Abbrev",,')
    for i in data:
        for j in data[i]:
            # j[0] = node1 type, j[1] = node1 name, j[2] = node2 type, j[3] = node2 name, j[4] = interaction sign (positive, negative, neutral), j[5] = node1 region, j[6] = node2 region
            writeMe.append('"general","'+i+'","'+j[0]+'","'+j[1]+'","'+j[2]+'","'+j[3]+'","'+j[4]+'","'+j[5]+'","'+j[6]+'",,')
    outFile = open(filename,'w')
    outFile.write('\n'.join(writeMe))
    outFile.close()

# Make an MDraw space delimited file
def mdraw(filename, data):
    node = 0
    nodeIds = {}
    writeMe1 = []
    for region in data:
        for j in data[region]:
            if not j[1] in nodeIds:
                node += 1
                nodeIds[j[1]] = node
            if not j[3] in nodeIds:
                node += 1
                nodeIds[j[3]] = node
            int1 = 0
            if j[4]=='positive':
                int1 = 1
            if j[4]=='negative':
                int1 = 2
            tmp1 = str(nodeIds[j[1]])+' '+str(nodeIds[j[3]])+' '+str(int1) #+' 1' # 
            if not tmp1 in writeMe1 and not j[1]==j[3]:
                writeMe1.append(tmp1)
    outFile = open(filename+'_mdraw.txt','w')
    outFile.write('\n'.join(writeMe1))
    outFile.close()
    outFile = open(filename+'_mdraw_LIST.txt','w')
    outFile.write('\n'.join([str(i)+' '+str(nodeIds[i]) for i in nodeIds]))
    outFile.close()


# Load TF->target gene dictionary
with open('tfbsDb_plus_and_minus_5000_entrez.json', 'r') as inFile:
    data = json.load(inFile)

# Load up id to motif thesaurus
id2Motif = {}
with open('id_conversion/humanTFs_All.csv','r') as inFile:
    header = inFile.readline().strip().split(',')
    while 1:
        inLine = inFile.readline()
        if not inLine:
            break
        split = inLine.strip().split(',')
        if not split[2] in id2Motif:
            id2Motif[split[2]] = []
        id2Motif[split[2]].append(split[0])

# To translate gene ids later
symbol2entrez = {}
with open('id_conversion/gene2entrezId.csv','r') as inFile:
    while 1:
        inLine = inFile.readline()
        if not inLine:
            break
        split = inLine.strip().split(',')
        symbol2entrez[split[1]] = split[0]

# TF family expansion
family2Id = {}
id2Family = {}
with open('id_conversion/tfFamilies.csv','r') as inFile:
    header = inFile.readline()
    while 1:
        inLine = inFile.readline()
        if not inLine:
            break
        split = inLine.split(',')
        split[2] = split[2].replace(' ',',').strip().split(',')
        family2Id[split[0]] = split[2]
        for splitId in split[2]:
            id2Family[splitId] = split[0]

            
#grab input TFs 
inputTfs = []
for tumortype in cluster2tf.keys():
    inputTfs = cluster2tf[tumortype]
    if tumortype is 'GBM':
        inputTfs = ['430', '1052', '1053', '1385', '84699', '9586', '1871', '1874', '144455', '79733', '1960', '1997', '2002', '2004', '80712', '2114', '2115', '2120', '51513', '2551', '2623', '2624', '2625', '9421', '3232', '10320', '3659', '3662', '3670', '91464', '3726', '10661', '11278', '128209', '10365', '9314', '1316', '51176', '9935', '23269', '4602', '4774', '4790', '7025', '9480', '5468', '5914', '5916', '3516', '5971', '864', '6257', '4093', '6659', '6660', '6662', '25803', '347853', '30009', '9496', '6929', '6925', '8463', '7022', '29842', '10155', '6935', '132625', '23051', '85416', '7707', '7764', '23528', '201516']
    # Load up exprssion data
    df = pd.read_csv('exprs_all/'+tumortype+'_RNAseq.csv', header=0, index_col=0)

    # Create a TFreg -> TFtarg dictionary
    # TFreg -> TFtarg
    # Type I Family Expansion = expand to all possible TFtargs via all motifs from family members
    # Type II Family Expansion = If no motif for TFreg then expand to family members
    print '  Building network: '+tumortype
    tfCascade = {}
    corTfCascade = {}
    output = {'CHIR_curve':[]}
    famExpType = args.family
    for TFreg in inputTfs:
        motifs = []
        if TFreg in id2Motif:
            motifs += id2Motif[TFreg]
        if TFreg in id2Family:
            if famExpType==1 or (famExpType==2 and len(motifs)==0):
                for TFregExp in family2Id[id2Family[TFreg]]:
                    if (not TFregExp==TFreg) and TFregExp in id2Motif:
                        motifs += id2Motif[TFregExp]
        # Iterate through motifs
        for motif in motifs:
            if motif in data:
                for geneTarg in data[motif]:
                    if geneTarg in inputTfs:
                        if not TFreg in tfCascade:
                            tfCascade[TFreg] = []
                        if not geneTarg in tfCascade[TFreg]:
                            tfCascade[TFreg].append(geneTarg)
                            if int(TFreg) in list(df.index) and int(geneTarg) in list(df.index): 
                                #new stuff to filter empties
                                TFregList = []
                                geneTargList = []        
                                for sizeIndex in range(df.shape[1]):
                                    check1 = df.loc[int(TFreg)][sizeIndex]
                                    check2 = df.loc[int(geneTarg)][sizeIndex]
                                    if not pd.isnull(check1):
                                        if not pd.isnull(check2):
                                            TFregList.append(check1)
                                            geneTargList.append(check2)
                                r1 = pearsonr(TFregList,geneTargList)
                                if abs(r1[0])>args.correlation and r1[1]<=0.05:
                                    print 'in',TFreg, geneTarg, r1
                                    if not TFreg in corTfCascade:
                                        corTfCascade[TFreg] = []
                                    if not geneTarg in corTfCascade[TFreg]:
                                        corTfCascade[TFreg].append(geneTarg)
                                    if r1[0]>0 and TFreg in symbol2entrez and geneTarg in symbol2entrez:
                                        output['CHIR_curve'].append(['gene',symbol2entrez[TFreg],'gene',symbol2entrez[geneTarg],'positive','A','A'])
                                    if r1[0]<0 and TFreg in symbol2entrez and geneTarg in symbol2entrez:
                                        output['CHIR_curve'].append(['gene',symbol2entrez[TFreg],'gene',symbol2entrez[geneTarg],'negative','A','A'])
    
    
    
    def down(key, tfCascade):
        """Function to identify all downstream targets of TF
        key.
    
        Args:
            key: starting node label.
            tfCascade: the full gene network as dict.
    
        Returns:
            A list of TF targets of key.
    
        """
        if key in tfCascade:
            print 'down', tfCascade[key]
            return tfCascade[key]
        else:
            return []
    
    def up(key,tfCascade):
        """Function to identify all upstream regulators of TF
        key.
    
        Args:
            key: starting node label.
            tfCascade: the full gene network as dict.
    
        Returns:
            A list of TF regulators of key.
    
        """
        outList = []
        for TFreg in tfCascade:
            if key in tfCascade[TFreg]:
                outList.append(TFreg)
        print 'up', outList
        return outList
    
    def subNetwork2(key,hops,tfCascade):
        """Function to grab out subnetwork from tfCascade given 
        a specific starting node (key) and for a given number of 
        hops.
    
        Args:
            key: starting node label.
            hops: number of node jumps away from starting node.
            tfCascade: the full gene network as dict.
    
        Returns:
            A dict of the subnetwork where the keys are TFregs
            and values are TFtargs.
    
        """
        mainList = [key]
        count = 0
        while count < hops:
            for g in mainList:
                tempList = []
                tempList += down(g,tfCascade)
                tempList += up(g,tfCascade)
                temp2List = []
                for x in tempList:
                    if not x in mainList:
                        temp2List.append(x)
            mainList = temp2List
            count += 1
        print mainList
    
        subNetwork = {}
        for TFreg in mainList:
            if not TFreg in subNetwork:
                subNetwork[TFreg] = []
            if TFreg in tfCascade:
                for TFtarg in tfCascade[TFreg]:
                    if TFtarg in mainList:
                        subNetwork[TFreg].append(TFtarg)
        return subNetwork
    
    #subGeneNetwork = subNetwork2(sPoint,hop,tfCascade)
    
    biotapestry('PanOutput/biotapestry_CHIR_curve'+tumortype+'.csv', output, ['CHIR_curve'])
    
    mdraw('PanOutput/mdraw_'+tumortype, output)

