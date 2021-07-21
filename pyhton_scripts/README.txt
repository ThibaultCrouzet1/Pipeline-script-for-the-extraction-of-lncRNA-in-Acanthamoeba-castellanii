The following scripts are used in the pipeline and were developed by me. The code of those script can be found on the following GitHub (link to the GitHub). Create the github


scriptGetFastaSeqById.py

This script uses the environment biopython and has 3 options: -t, -i, -o
The option -t corresponds to the fasta file of the transcriptome.
The option -i corresponds to id list fille.
The option -o corresponds to the output file.

The goal of this script is to extract the sequences of a list of IDs into a fasta file.


merge_annot.py

This script uses the environment jupyterPandas and has 3 merging functions chosen with the option --merging : 
1.	merging_table with –merging = "without_model"
2.	merging_table_output_with_PLEK_model_output –merging = "with_model"
3.	merging_table_output_with_PLEK_model_balanced_output –merging = "with_model_balanced"

The first function, merging_table, takes the annotation files from blastx and InterPro to merge them with the prediction result from CPC2 and PLEK without any model.
The second function, merging_table_ouput_with_PLEK_model_output, combines the result of the first function after a formating of the columns (space to tab) with the prediction results from PLEK with the model. 
The third function, merging_table_output_with_PLEK_model_balanced_output, combines the result of the second function after a formating of the columns (space to tab) with the result of the PLEK prediction with a balanced model.

The options of the function are:
•	--merging, this option allows to choose which function of the script is used.
•	--tr_inputfile, this option represents the labelled transcriptome in the fasta file format  
•	--mergedtab, this option represents the merged table created by the first and second function
•	--blastx, this option represents the list file of the id of the non-annotated transcripts from the blastx annotation file
•	--interpro, this option represents the condensed list file of the interPro scan annotation where the condensed list takes the id, the domain and prediction done by InterPro scan
•	--cpc2, this option represents the prediction file of the lncRNA done by CPC2
•	--plek, this option represents the prediction file of the lncRNA done by PLEK without any model
•	--plekmodel, this option represents the prediction file of the lncRNA done by PLEK with a model (balanced or not)
•	--output, this option represents the name of the output file

The first function takes the options:
--tr_inputfile, --blastx, --interpro, --cpc2, --plek, --output

The second function takes the options:
--mergedtab, --plekmodel, --output
Where the --plekmodel option take the PLEK prediction with the unbalanced model.

The third function takes the options:
--mergedtab, --plekmodel, --output
Where the --plekmodel option take the PLEK prediction with the balanced model.



extraction.py

This script uses the environment jupyterPandas and its goal is to extract the lncRNA from the merged table produced by the third function of the merged_annot.py script into a data frame with the id, the length of the transcript, the annotations, the prediction data.
A transcript is considered a lncRNA if:
•	the length of the transcript is 200 nucleotides or more, 
•	both annotations (blastx and interPro) describe the transcript noncoding (no hit) 
•	the prediction of CPC2 is non-coding and the CPC2 prediction score is less than 0,4 
•	the PLEK predictions with the balanced model and without any model are non-coding with both a prediction score inferior to 0 

The script also extracts the sequence of the lncRNA extracted and add it to their transcript line inside the data frame.  

The script has 3 options: 
•	-m, the merged table obtain with the third function of the merge_annot.py script
•	-t, the labelled transcriptome in the fasta file format to extract the sequence of the lncRNA
•	-o, the name of the output file.



getExpressionLevel.py

This script extracts the expression level of the lncRNA from the expression level files of Ac and AcPa. And uses the biopython environment.
The script takes as input the expression level file and the fasta sequences of the extracted lncRNA.
The script has 3 options:
•	-f, the lncRNA fasta file
•	-e, the expression level file
•	-o, output file name



sRNA_extraction.py

This script uses the environment jupyterPandas and takes the merged table and the labelled transcriptome in the fasta file format to extract the sequence of the lncRNA as an input.
The script extracts all the transcripts with their annotations and sequence that have a length inferior to 200 nucleotides.

The script has 3 options: 
•	-m, the merged table obtain with the third function of the merge_annot.py script
•	-t, the labelled transcriptome in the fasta file format to extract the sequence of the lncRNA
•	-o, the name of the output file.


sRNA_merge.py

This script uses the environment jupyterPandas and its goal is to merge all the prediction files and annotations.
The options of the function are:
•	--sRNA_inputfile, this option represents the labelled transcriptome in the fasta file format  
•	--cpc2, this option represents the prediction file of the lncRNA done by CPC2
•	--plek, this option represents the prediction file of the lncRNA done by PLEK without any model
•	--plekmodel, this option represents the prediction file of the lncRNA done by PLEK with a model non-balanced
•	--plekmodel_balanced, this option represents the prediction file of the lncRNA done by PLEK with a model balanced
•	--output, this option represents the name of the output file



seqInBothTr.py

This script uses the environment biopython and its goal is to compare, and extract expressed lncRNA sequences that are in both Ac and AcPa.
The script takes as input the fasta file of the expressed lncRNA from Ac and AcPa.
This script has the following option:
•	-a, the fasta file of the Ac expressed lncRNA
•	-p, the fasta file of the AcPa expressed lncRNA
•	-o, the name of the output file
