PanCancer-TFNetwork-Motifs
Using Python transcription factor networks are generated across multiple cancer types. Command line version of FANMOD pulls enriched motifs which are analysed and plotted with survival data.

Wrapper file can be run via command line given a correlation input value:
python PanCancerTFNetworkMotifWrapper.py --correlation 0.5
Input files are limited to only running 'Example' tumor type. This can be expanded by modifying the input list of tumor types in the wrapper file.
Default TF family expansion is ON. This can be changed by increasing the integer in the wrapper file from 1 to 2.


File 1)
regulators_vSurv.py
Used to generate a set number of input genes for gene network generation.
Input Files:
1)gbm/gbm_includedBiclusters.csv
-list of bioclusters
2)id_conversion/gene2entrez.csv
-Translation
3)gbm/mature_hsa
-ID to accession number
4)postProcessed_vSurv/postProcessed_tumorType_pita.csv
-Biocluster, Genes, Patients, Tumor Expression, various Tumor Attributes
OutputFiles:
1)RegulatorsOutput/output_vSurv.csv
-Bioclusters from TFs


File 2)
TfCascade.py
Generates the TF Network
Input Files:
1)RegulatorsOutput/output_vSurv.csv
-From File 1, Input TFs for Network
2)tfbsDb_plus_and_minus_5000_entrez.json
-Data to build of networks
3)id_conversion/humanTFs_All.csv
-Motif Name, New Symbol, Entrex ID
4)id_conversion/gene2entrezId.csv
-Translation
5)id_conversion/tfFamilies.csv
-Grouped Genes to a TF family
Output Files:
1)PanOutput/biotapestry_CHIR_curvetumortype.csv
-Source to target network information
2)PanOutput/mdraw_Tumor_mdraw.txt
-Source to Target, Sign (for FANMOD)
3)PanOutput/mdraw_Tumor_mdraw_LIST.txt
-Converts mdraw gene numbers to gene symbols

File 3)
fanmod_command_line_linux
Finds enriched motifs in networks
Input Files:
1)PanOutput/mdraw_Tumor_mdraw.txt
-From File 2, Source to target nodes with sign
Output Files:
1)fanmod_motifs_Tumor.txt
-Output of enriched motifs
2)fanmod_motifs_Tumor.txt.dump
-adjency matrix with participating vertices, easier to parse
For more information and files: https://github.com/gtremper/Network-Motif/tree/master/fanmod/FANMOD-command_line-source

File 4)
consistilatorV2.py
Finds correlation between nodes of each motif.
Input Files:
1)exprs_all/Tumor_RNAseq.csv
-Tumor gene expression data
2)FanmodOutput/fanmod_motifs_Tumor.txt
-From File 3, output of enriched motifs
3)PanOutput/mdraw_Tumor_mdraw_LIST.txt
-From File 2, converts mdraw numbers to gene symbols
4)FanmodOutput/fanmod_motifs_Tumor.txt.dump
-From File 3, adjency matrix with participating vertices
5)PanOutput/biotapestry_CHIR_curveTumor.csv
-From File 2, Source to target network information
6)id_conversion/gene2entrezId.csv
-Translation
Output Files:
1)ConsistanceResults/results_Subset_Tumor.csv
-Motifs with consistency information

File 5)
plotNetworkMotifs.py
Plots information from previous files to PDF
Input Files:
1)PanOutput/biotapestry_CHIR_curveTumors.csv
-From File 2, Source to target network information
2)ConsistanceResults/results_Subset_Tumor.csv
-Motifs with consistency information
3)id_conversion/gene2entrezId.csv
-Translation
4)exprs_all/Tumor_RNAseq.csv
-Tumor gene expression data
5)phenotypes_panCancer.csv
-Patient measurements
Output Files:
1)PipeOut/TumorsNetMotifs_3node_seq.pdf
-pfd output of figs
2)PipeOut/pairwise_comparisons.csv
-Motifs full information
