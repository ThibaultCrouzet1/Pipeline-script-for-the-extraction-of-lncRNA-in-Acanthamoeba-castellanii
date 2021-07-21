import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import glob, os, sys, getopt

#function
def getFastaFromAnnot(merged_table,seq_fasta,output):
    #openning the megre table in a dataframe
    print("Creation of the dataframe for analysis ...")
    df_input = pd.read_table(merged_table,dtype={'Seq_Length': str,'CPC2_score' : str},names=['ligne','ID','Seq_Length','Annotation_Balstx','Annotation_Interpro_TSV','CPC2_score','CPC2_prediction','PLEK_score','PLEK_prediction','PLEK_score_model','PLEK_prediction_model','PLEK_score_model_balanced','PLEK_prediction_model_balanced'] )
    df_input = df_input.drop(columns="ligne")
    df_input = df_input.iloc[1:]
    print("Complete ...")

    #opening the transcript file seq_fasta
    print("Getting the sequences ...")
    records = list(SeqIO.parse( seq_fasta, "fasta")) #validated

    #creating the output fasta file for coding and non-coding blastx and interproscan 
    ncBlastx = open(output+"blastx_Non_coding.fasta",'w')
    cBlastx = open(output+"blastx_Coding.fasta",'w')
    ncInterpro = open(output+"interproscan_Non_coding.fasta",'w')
    cInterpro = open(output+"interproscan_Coding.fasta",'w')


    #for each row of the merge table dataframe
    print("Analysing and extracting the fasta ...")
    for index, row in df_input.iterrows() :
        #float() to convert the str value to their scientific numerical value (float)
        #if both annotation file say non coding​ + CPC score < 0,4​ + PLEK score < 0​ or PLEK score with the model < 0 + the prediction of both PLECK or PLEK with the model & CPC2 = non coding​ +  the length of the seq >= 200nt
        #for the id of the seq we get the seq from the seq_fasta file with the SeqIO module when the id of the seq equals the one in the row 
        for record in records :
            if record.id == row.loc['ID'] :
                sequence = str(record.seq)
                break

        if (row.loc['Annotation_Balstx'] == "Non_coding" ) :
            ncBlastx.write(">"+row.loc['ID']+"\n"+sequence+"\n")

        elif (row.loc['Annotation_Balstx'] == "Coding" ) :
            cBlastx.write(">"+row.loc['ID']+"\n"+sequence+"\n")

        if (row.loc['Annotation_Interpro_TSV'] == "non_coding") or (row.loc['Annotation_Interpro_TSV'] == "Non_Coding") :
            ncInterpro.write(">"+row.loc['ID']+"\n"+sequence+"\n")

        elif (row.loc['Annotation_Interpro_TSV'] == "coding") :
            cInterpro.write(">"+row.loc['ID']+"\n"+sequence+"\n")
        

    ncBlastx.close()
    cBlastx.close()
    ncInterpro.close()
    cInterpro.close()

    print("Output file created.")





#calling the function and its parameters
def main(argv):
    input_outputmergetabfile = ''
    tr_label_file = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"m:t:o:",["mergedtab=","tr_inputfile=","outputfiles="])

    except getopt.GetoptError:
        print('getFastaTrinityStats.py -m <output_merge_model_format.tab> -t <tr_label.fasta> -o <outputfile>')
        print('getFastaTrinityStats.py --mergedtab=<output_merge_model_format.tab> --tr_inputfile=<tr_label.fasta> --outputfile=<outputfiles>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('getFastaTrinityStats.py -m <output_merge_model_format.tab> -t <tr_label.fasta> -o <outputfile>')
            print('getFastaTrinityStats.py --mergedtab <output_merge_model_format.tab> --tr_inputfile <tr_label.fasta> --outputfile <outputfile>')
            sys.exit()
        elif opt in ("-m", "--mergedtab"):
            input_outputmergetabfile = arg
        
        elif opt in ("-t", "--tr_inputfile"):
            tr_label_file = arg
        elif opt in ("-o", "--outputfile"):
            outputfile = arg
			

    getFastaFromAnnot(input_outputmergetabfile,tr_label_file,outputfile)

if __name__ == "__main__":
   main(sys.argv[1:])

#example of parameters
#getFastaTrinityStats.py -m "04.LNCRNA_PREDICTION/4.3.results_observations/with_model/Ac/Ac_output_merge_model_balanced_format.tab" -t "01.LABELED_DATA/Ac_labeled/Ac_tr_label.fasta" -o "04.LNCRNA_PREDICTION/4.3.results_observations/Ac_"
#getFastaTrinityStats.py -m "04.LNCRNA_PREDICTION/4.3.results_observations/with_model/AcPa/AcPa_output_merge_model_balanced_format.tab" -t "01.LABELED_DATA/AcPa_labeled/AcPa_tr_label.fasta" -o "04.LNCRNA_PREDICTION/4.3.results_observations/AcPa_"


