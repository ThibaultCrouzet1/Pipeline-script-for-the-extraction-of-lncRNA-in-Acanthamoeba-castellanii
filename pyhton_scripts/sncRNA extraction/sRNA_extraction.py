import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import glob, os, sys, getopt


def extractSRNA(merged_table,seq_fasta,output):
    #openning the megre table in a dataframe
    print("Creating the dataframe for analysis ...")
    df_input = pd.read_table(merged_table,dtype={'Seq_Length': str,'CPC2_score' : str},names=['ligne','ID','Seq_Length','Annotation_Blastx','Annotation_Interpro_TSV','CPC2_score','CPC2_prediction','PLEK_score','PLEK_prediction','PLEK_score_model','PLEK_prediction_model','PLEK_score_model_balanced','PLEK_prediction_model_balanced'] )
    df_input = df_input.drop(columns="ligne")
    df_input = df_input.iloc[1:]
    print("Complete ...")
    
    #opening the transcript file seq_fasta
    print("Getting the sequences ...")
    records = list(SeqIO.parse( seq_fasta, "fasta")) 
    
    
    #creation of a new empty dataframe identical to the merged table one but with a column for the sequence 
    print("Preparation for the output file")
    df_output = pd.DataFrame(columns=['ID','Seq_Length','Annotation_Blastx','Annotation_Interpro_TSV','Sequence'])
    list_columns=['ID','Seq_Length','Annotation_Blastx','Annotation_Interpro_TSV','Sequence']
    

    #for each row of the merge table dataframe
    print("Analysing and extracting the lncRNA ...")
    for index, row in df_input.iterrows() :
        #float() to convert the str value to their scientific numerical value (float)
        #if both annotation file say non coding​ + CPC score < 0,4​ + PLEK score < 0​ or PLEK score with the model < 0 + the prediction of both PLECK or PLEK with the model & CPC2 = non coding​ +  the length of the seq >= 200nt
        if ( int(row.loc['Seq_Length'])  < 200  ) :    
            #for the id of the seq we get the seq from the seq_fasta file with the SeqIO module and insert it into the row 
                record = next((r for r in records if r.id == row.loc['ID']))
                sequence = str(record.seq)

                #copy the row into a new row of the empty data frame
                output_row = pd.DataFrame( [[row.loc['ID'] , row.loc['Seq_Length'] , row.loc['Annotation_Blastx'] , row.loc['Annotation_Interpro_TSV'] , sequence  ]], columns=list_columns )
                df_output = df_output.append( output_row , ignore_index = True ) # ignore_index = True allow to have only one index and not to fuse both index
    print("Complete ...")
    #writing the output dataframe inside a file
    #opening the output
    out=open(output+".tab","w")
    out.write(df_output.to_string())

    print("Extracting the sequences into fasta file ...")
    outfasta = open(output+".fasta","w")
    records = list(SeqIO.parse( seq_fasta, "fasta")) 
    for r200 in records :
        if( len(r200.seq) < 200 ):
            sequence = str(r200.seq)
            print(sequence)
            outfasta.write(">"+r200.id+" length "+str(len(r200.seq))+"\n"+sequence+"\n")
    

    print("Output files created.")
    print(df_output.head())





#calling the function and its parameters
def main(argv):
	input_outputmergetabfile = ''
	tr_label_file = ''
	outputfile = ''
	try:
		opts, args = getopt.getopt(argv,"m:t:o:",["mergedtab=","tr_inputfile=","outputfile="])
		
	except getopt.GetoptError:
		print('sRNA_extraction.py -m <output_merge_model_format.tab> -t <tr_label.fasta> -o <outputfile.tab>')
		print('sRNA_extraction.py --mergedtab=<output_merge_model_format.tab> --tr_inputfile=<tr_label.fasta> --outputfile=<outputfile.tab>')
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print('sRNA_extraction.py -m <output_merge_model_format.tab> -t <tr_label.fasta> -o <outputfile.tab>')
			print('sRNA_extraction.py --mergedtab <output_merge_model_format.tab> --tr_inputfile <tr_label.fasta> --outputfile <outputfile.tab>')
			sys.exit()
		elif opt in ("-m", "--mergedtab"):
			input_outputmergetabfile = arg
			
		elif opt in ("-t", "--tr_inputfile"):
			tr_label_file = arg
			
		elif opt in ("-o", "--outputfile"):
			outputfile = arg
			

	extractSRNA(input_outputmergetabfile,tr_label_file,outputfile)

if __name__ == "__main__":
   main(sys.argv[1:])


#example of parameters
#-m 04.LNCRNA_PREDICTION/4.3.results_observations/with_model/Ac/Ac_output_merge_model_balanced_format.tab -t 01.LABELED_DATA/Ac_labeled/Ac_tr_label.fasta -o 07.SRNA_EXTRACTION/Ac/Ac_sRNA_extracted

#-m 04.LNCRNA_PREDICTION/4.3.results_observations/with_model/AcPa/AcPa_output_merge_model_balanced_format.tab -t 01.LABELED_DATA/AcPa_labeled/AcPa_tr_label.fasta -o 07.SRNA_EXTRACTION/AcPa/AcPa_sRNA_extracted
