import subprocess
import argparse

parser = argparse.ArgumentParser(description='use for testing argparse.py')
parser.add_argument('--correlation', help='Run at a specific correlation', type = float)
args = parser.parse_args()

#list of tumor types to run. Trial set to ACC.
tumors = 'ACC'#,BLCA,BRCA,CESC,CHOL,COAD,DLBC,ESCA,GBM,HNSC,KICH,KIRC,KIRP,LGG,LIHC,LUAD,LUSC,MESO,OV,PAAD,PCPG,PRAD,READ,SARC,SKCM,STAD,TGCT,THCA,THYM,UCEC,UCS,UVM'

correlation = args.correlation
#correlation = 0.5

# Type I Family Expansion = expand to all possible TFtargs via all motifs from family members
# Type II Family Expansion = If no motif for TFreg then expand to family members
family = 1


#Communicate to Regulators-vSurv
a = 'python regulators_vSurv.py  --correlation '+str(correlation)+' --tumors '+str(tumors)
sOne = subprocess.Popen(a, shell=True)
sOne.communicate()
    
#Communicate to TfCascade
b = 'python TfCascade.py  --correlation '+str(correlation)+' --family '+str(family)
sTwo = subprocess.Popen(b, shell=True)
sTwo.communicate()

#Communicate to FANMOD
tumors = [item for item in tumors.split(',')]
for tumorType in tumors:
    c = './fanmod_command_line_linux 3 100000 1 PanOutput/mdraw_'+str(tumorType)+'_mdraw.txt 1 0 1 2 0 1 0 1000 3 3 FanmodOutput/fanmod_motifs_'+str(tumorType)+'.txt 1 1'
    sThree = subprocess.Popen(c, shell=True)
    sThree.communicate()

#communicate to ConsistilatorV2
tumors = [item for item in tumors.split(',')]
for tumorType in tumors:
    d = 'python consistilatorV2.py  --tumor '+str(tumorType)
    sThree = subprocess.Popen(d, shell=True)
    sThree.communicate()

#communicate to plotMotifs
tumors = [item for item in tumors.split(',')]
for tumorType in tumors:
    e = 'python plotNetworkMotif.py  --tumors '+str(tumorType)
    sThree = subprocess.Popen(e, shell=True)
    sThree.communicate()