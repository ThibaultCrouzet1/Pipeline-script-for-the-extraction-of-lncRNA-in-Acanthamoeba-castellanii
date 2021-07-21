import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import glob, os, sys, getopt


def extractLNCRNA(merged_table,seq_fasta,output):
    #openning the megre table in a dataframe
    print("Creation of the dataframe for analysis ...")
    df_input = pd.read_table(merged_table,dtype={'Seq_Length': str,'CPC2_score' : str},names=['ligne','ID','Seq_Length','Annotation_Balstx','Annotation_Interpro_TSV','CPC2_score','CPC2_prediction','PLEK_score','PLEK_prediction','PLEK_score_model','PLEK_prediction_model','PLEK_score_model_balanced','PLEK_prediction_model_balanced'] )
    df_input = df_input.drop(columns="ligne")
    df_input = df_input.iloc[1:]
    print("Complete ...")
    
    #opening the transcript file seq_fasta
    print("Getting the sequences ...")
    records = list(SeqIO.parse( seq_fasta, "fasta")) #validated
    
    
    #creation of a new empty dataframe identical to the merged table one but with a column for the sequence 
    print("Preparation for the output file")
    df_output = pd.DataFrame(columns=['ID','Seq_Length','Annotation_Balstx','Annotation_Interpro_TSV','CPC2_score','CPC2_prediction','PLEK_score','PLEK_prediction','PLEK_score_model','PLEK_prediction_model','PLEK_score_model_balanced','PLEK_prediction_model_balanced','Sequence'])
    list_columns=['ID','Seq_Length','Annotation_Balstx','Annotation_Interpro_TSV','CPC2_score','CPC2_prediction','PLEK_score','PLEK_prediction','PLEK_score_model','PLEK_prediction_model','PLEK_score_model_balanced','PLEK_prediction_model_balanced','Sequence']
    
    #for each row of the merge table dataframe
    print("Analysing and extracting the lncRNA ...")
    for index, row in df_input.iterrows() :
        #float() to convert the str value to their scientific numerical value (float)
        #if both annotation file say non coding​ + CPC score < 0,4​ + PLEK score < 0​ or PLEK score with the model < 0 + the prediction of both PLECK or PLEK with the model & CPC2 = non coding​ +  the length of the seq >= 200nt
        if ( int(row.loc['Seq_Length'])  >= 200  ) and ( (row.loc['Annotation_Balstx'] == "Non_coding" ) and (row.loc['Annotation_Interpro_TSV'] == "non_coding") or (row.loc['Annotation_Interpro_TSV'] == "Non_Coding") ) :
            if (( float(row.loc['CPC2_score']) < 0.40 ) and ( row.loc['CPC2_prediction'] == "noncoding" )) and ( ((float(row.loc['PLEK_score']) < 0 ) and (row.loc['PLEK_prediction'] == "Non-coding")) and ((float(row.loc['PLEK_score_model_balanced']) < 0 ) and (row.loc['PLEK_prediction_model_balanced'] == "Non-coding")) ) :
                
                #for the id of the seq we get the seq from the seq_fasta file with the SeqIO module and insert it into the row 
                record = next((r for r in records if r.id == row.loc['ID']))
                sequence = str(record.seq)

                #copy the row into a new row of the empty data frame
                output_row = pd.DataFrame( [[row.loc['ID'] , row.loc['Seq_Length'] , row.loc['Annotation_Balstx'] , row.loc['Annotation_Interpro_TSV'] , row.loc['CPC2_score'] , row.loc['CPC2_prediction'] , row.loc['PLEK_score'] , row.loc['PLEK_prediction'] , row.loc['PLEK_score_model'] , row.loc['PLEK_prediction_model'], row.loc['PLEK_score_model_balanced'] , row.loc['PLEK_prediction_model_balanced'],sequence  ]], columns=list_columns )
                df_output = df_output.append( output_row , ignore_index = True ) # ignore_index = True allow to have only one index and not to fuse both index
    print("Complete ...")
    #writing the output dataframe inside a file
    #opening the output
    out=open(output,"w")
    out.write(df_output.to_string())
    print("Output file created.")
    print(df_output.head())


#calling the function and its parameters
def main(argv):
    input_outputmergetabfile = ''
    tr_label_file = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"m:t:o:",["mergedtab=","tr_inputfile=","outputfile="])

    except getopt.GetoptError:
        print('extraction.py -m <output_merge_model_format.tab> -t <tr_label.fasta> -o <outputfile.tab>')
        print('extraction.py --mergedtab=<output_merge_model_format.tab> --tr_inputfile=<tr_label.fasta> --outputfile=<outputfile.tab>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('extraction.py -m <output_merge_model_format.tab> -t <tr_label.fasta> -o <outputfile.tab>')
            print('extraction.py --mergedtab <output_merge_model_format.tab> --tr_inputfile <tr_label.fasta> --outputfile <outputfile.tab>')
            sys.exit()
        elif opt in ("-m", "--mergedtab"):
            input_outputmergetabfile = arg
        
        elif opt in ("-t", "--tr_inputfile"):
            tr_label_file = arg
        elif opt in ("-o", "--outputfile"):
            outputfile = arg
			

    extractLNCRNA(input_outputmergetabfile,tr_label_file,outputfile)

if __name__ == "__main__":
   main(sys.argv[1:])



#example of parameters
#extraction.py -m "04.LNCRNA_PREDICTION/4.3.results_observations/with_model/Ac/Ac_output_merge_model_balanced_format.tab" -t "01.LABELED_DATA/Ac_labeled/Ac_tr_label.fasta" -o "05.LNCRNA_EXTRACTION/5.1.EXTRACTED_LNCRNA/Ac/Ac_lncRNA_extracted_balanced.tab"
#extraction.py -m "04.LNCRNA_PREDICTION/4.3.results_observations/with_model/AcPa/AcPa_output_merge_model_balanced_format.tab" -t "01.LABELED_DATA/AcPa_labeled/AcPa_tr_label.fasta" -o "05.LNCRNA_EXTRACTION/5.1.EXTRACTED_LNCRNA/AcPa/AcPa_lncRNA_extracted_balanced.tab"
