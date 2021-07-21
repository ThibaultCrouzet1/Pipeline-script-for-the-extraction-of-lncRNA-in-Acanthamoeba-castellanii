import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import glob, os, sys, getopt

def merging_table(transcript, blastx_file,interpro_TSV_file,filecpc2, fileplek,output):
    #opening the transcript file
    records = list(SeqIO.parse( transcript, "fasta"))
    #opening the output
    out=open(output,"w")
    
    #base dataframe for the merging that contians the id and length of the seq 
    print("Creation of the base dataframe for the merging that contains the id and length for all the transcript sequences ...")
    seq_tr_list=list()
    for record in records:
        #get the id + the length of the seq
        seq_tr_list.append([record.id,len(record.seq)])
    

    #creation of the dataframe with only the id of the seq and its length
    df = pd.DataFrame(seq_tr_list,columns=['ID','Seq_Length'])
    print("Dataframe created.\nAnalysis of the annotation files starting ...")


    #creation of the dataframe for the annotation from Blastx
    

    #we have a file that contains the non coding sequence id from blastx file 
    #so if we take the id from that are in the file we can say that in the column annot blastx the corresponding id is noncoding if present or coding if absent
    print("Creation of the Blastx annotation dataframe ...")
    annot_blastx_file = open(blastx_file, 'r')
    blastx_id_list=list()
    for line_blastx in annot_blastx_file.readlines():
        id_blastx = line_blastx.split(' ')[0]
        blastx_id_list.append(id_blastx)
    df_blastx=pd.DataFrame(blastx_id_list,columns={'ID'})
    df_blastx=df_blastx.assign(Annotation_Balstx='Non_coding')
    print("Dataframe created ...")
    



    #creation of the dataframe for the annotation from Interprot
    #for the TSV (tab) file
    #command bash that will extract all the id of the sequence that have MobiDBLite and consensus disorder prediction on the line and then put the reuslt into a file

    #create a dataframe from the file with a column annoted_interpro where the default value is non coding or non coding
    #preparation of the data into a non filtered dataframe
    print("Preparation for the creation of the interpro TSV annotation dataframe ...")
    annot_interpro_TSV_file = open(interpro_TSV_file, 'r')
    interpro_TSV_id_list=list()
    for line_interpro_TSV in annot_interpro_TSV_file.readlines():
        id_interpro_TSV = line_interpro_TSV.split('\t')
        if id_interpro_TSV[3] == '\n' :
            interpro_pred = [id_interpro_TSV[0],id_interpro_TSV[2],"\n"] #the \n will serve as a beacon to avoid empty row inside the column
        else :
            id_interpro_TSV[3] = id_interpro_TSV[3].replace(' ', '_')
            id_interpro_TSV[3] = id_interpro_TSV[3].replace('\n', '')
            interpro_pred = [id_interpro_TSV[0],id_interpro_TSV[2],id_interpro_TSV[3]]
        
        interpro_TSV_id_list.append(interpro_pred)
    
    df_interpro_TSV_non_filtered=pd.DataFrame(interpro_TSV_id_list, columns=['ID', 'DOMAIN', 'PREDICTION'])
    
    #get all the different id only once to avoid duplicate
    ID_list = df_interpro_TSV_non_filtered['ID'].unique()

    #setting the index to speed up the loc function
    df_interpro_TSV_non_filtered=df_interpro_TSV_non_filtered.set_index(['ID'])
    

    #filtration step
    print("Extracting the annotation from the TSV interpro output file ...")
    
    #creation of the filtered dataframe for the TSV interpro files
    coding_id_list = list()
    
    for id in ID_list :
        verif = False #confirm that for the id there is no coding prediction
        df_all_info_nofil = df_interpro_TSV_non_filtered.loc[[id]]
        for i in range(len(df_all_info_nofil)) :
            if df_all_info_nofil['DOMAIN'].iloc[i] != 'mobidb-lite' and  df_all_info_nofil['PREDICTION'].iloc[i] != 'consensus_disorder_prediction' :
                coding_id_list.append([id,'coding'])
                verif=True #confirm that for the id there is at least one coding prediction
                break
        if verif == False : #verify the confirmation that for the id there is no coding prediction
            coding_id_list.append([id,'non_coding'])
       
    print("Creation of the interpro TSV annotation dataframe ...")
    df_interpro_TSV=pd.DataFrame(coding_id_list,columns={'ID','Annotation_Interpro_TSV'})
    print(df_interpro_TSV.head())
    pd.set_option('display.max_rows',None)
    print("Dataframe created ...")
    
    
    
    print("Analysis of the prediction software output starting ...")
    print("CPC2 analysis starting.\nCreation of the dataframe commencing ...")
    #getting all the results of the cpc2 script into a dataframe
    df_cpc = pd.read_table(filecpc2)
    #droping the non essential columns
    df_cpc = df_cpc.drop(['transcript_length','peptide_length','Fickett_score','pI','ORF_integrity'], axis=1)
    #renaming the #ID columns to ID for the merging and to specify the origin
    df_cpc=df_cpc.rename(columns={"#ID":"ID","coding_probability":"CPC2_score","label":"CPC2_prediction"})
    print("CPC2 dataframe created ...")
    
    
    #getting all the results of the PLEK script into a dataframe
    #need to read the full file because of the structure of the id column
    print("Plek analysis starting.\nCreation of the dataframe commencing ...")
    plek = open(fileplek, 'r')
    pleklist=list()
    for lineplek in plek.readlines():
        fullidplek =  lineplek.split('\t')[2]
        idfastaplek =  fullidplek.split(' ')[0]
        idplek = idfastaplek.split('>')[1]
        pleklist.append([idplek,lineplek.split('\t')[1],lineplek.split('\t')[0]])
    df_plek= pd.DataFrame(pleklist,columns=['ID','PLEK_score','PLEK_prediction'])
    print("PLEK dataframe created ...")
    
    #merging all the dataframe into one
    print("Dataframe merging started ...")
    df_output = pd.merge(df,df_blastx, how='outer').fillna('Coding')
    df_output = pd.merge(df_output,df_interpro_TSV, left_on='ID', right_on='ID' ,how='outer') #after 10 000 seq the ID and the PRediction are inverted and cause the merging to failed
    df_output = pd.merge(df_output,df_cpc, how='outer')
    df_output = pd.merge(df_output,df_plek, how='outer')
    print("Merging finished ...")
    
    print("Preview of the merged dataframe ...")
    print(df_output)
    
    #writing of the output file 
    out.write(df_output.to_string())
    



def merging_table_ouput_with_PLEK_model_output(merged_table, fileplek, output):
    #openning the megre table in a dataframe
    print("Creation of the dataframe of the merged table ...")
    df_input = pd.read_table(merged_table,dtype={'Seq_Length': str,'CPC2_score' : str},names=['ligne','ID','Seq_Length','Annotation_Balstx','Annotation_Interpro_TSV','CPC2_score','CPC2_prediction','PLEK_score','PLEK_prediction'])
    df_input = df_input.drop(columns="ligne")
    df_input = df_input.iloc[1:]
    print("Creation completed ...")

    print("Plek with model analysis starting.\nCreation of the dataframe commencing ...")
    plek = open(fileplek, 'r')
    pleklist=list()
    for lineplek in plek.readlines():
        fullidplek =  lineplek.split('\t')[2]
        idfastaplek =  fullidplek.split(' ')[0]
        idplek = idfastaplek.split('>')[1]
        pleklist.append([idplek,lineplek.split('\t')[1],lineplek.split('\t')[0]])
    df_plek= pd.DataFrame(pleklist,columns=['ID','PLEK_score_model','PLEK_prediction_model'])
    print("PLEK dataframe created ...")

    #merging all the dataframe into one
    print("Dataframe merging started ...")
    df_output = pd.merge(df_input,df_plek, how='outer').fillna('Non_Coding')
    print("Merging finished ...")

    print("Preview of the merged dataframe ...")
    print(df_output)

    
    #opening the output
    out=open(output,"w")
    #writing of the output file 
    out.write(df_output.to_string())
    print("output file completed")


#take the merged table from the function "merging_table_ouput_with_PLEK_model_output" with the prediction result from PLEK with the balanced model
def merging_table_ouput_with_PLEK_model_balanced_output(merged_table, fileplek, output):
    #openning the megre table in a dataframe
    print("Creation of the dataframe of the merged table ...")
    df_input = pd.read_table(merged_table,dtype={'Seq_Length': str,'CPC2_score' : str},names=['ligne','ID','Seq_Length','Annotation_Balstx','Annotation_Interpro_TSV','CPC2_score','CPC2_prediction','PLEK_score','PLEK_prediction','PLEK_score_model','PLEK_prediction_model'])
    df_input = df_input.drop(columns="ligne")
    df_input = df_input.iloc[1:]
    print("Creation completed ...")

    print("Plek with model analysis starting.\nCreation of the dataframe commencing ...")
    plek = open(fileplek, 'r')
    pleklist=list()
    for lineplek in plek.readlines():
        fullidplek =  lineplek.split('\t')[2]
        idfastaplek =  fullidplek.split(' ')[0]
        idplek = idfastaplek.split('>')[1]
        pleklist.append([idplek,lineplek.split('\t')[1],lineplek.split('\t')[0]])
    df_plek= pd.DataFrame(pleklist,columns=['ID','PLEK_score_model_balanced','PLEK_prediction_model_balanced'])
    print("PLEK dataframe created ...")

    #merging all the dataframe into one
    print("Dataframe merging started ...")
    df_output = pd.merge(df_input,df_plek, how='outer').fillna('Non_Coding')
    print("Merging finished ...")

    
    print("Preview of the merged dataframe ...")
    print(df_output)

    
    #opening the output
    out=open(output,"w")
    #writing of the output file 
    out.write(df_output.to_string())
    print("output file completed")









#calling the function and its parameters
def main(argv):
    option=''
    mergedtab=''
    tr_inputfile=''
    blastx=''
    interpro=''
    cpc2=''
    plek=''
    plekmodel=''
    output=''
	
    try:
        opts, args = getopt.getopt(argv,"o",["merging=","mergedtab=","tr_inputfile=","blastx=", "interpro=", "cpc2=", "plek=", "plekmodel=", "output="])
		
    except getopt.GetoptError: 
        print("This script has 3 functions the first function take the annotation file from blastx and interpro combine them with the prediction result from CPC2 and PLEK(without the model).\nThe second function merging_table_ouput_with_PLEK_model_output combine the result of the last function after a formating of the columns (space to tab) with the prediction results from PLEK with the model. The third function take the result of second and add the result of the PLEK prediction with a balanced model.\nYou can choose the function to use with the option merging like followed :")
        print('merge_annot.py --merging="without_model" --tr_inputfile=<tr_label.fasta> --blastx=<non_annot_Ac_seqid.list> --interpro=<Ac_interpro_TSV.output.condense.list> --cpc2=<tr_label_cpc2_r.tab> --plek=<Ac_output_PLEK_tr4_ml200_k5.txt> --output=<outputfile.tab>')
        print('merge_annot.py --merging="with_model" --mergedtab=<Ac_output_merge_without_model_format.tab> --plekmodel=<Ac_output_PLEK_model_Dic_tr4_ml200_k5.txt> --output=<outputfile.tab>')
        print('merge_annot.py --merging="with_model_balanced" --mergedtab=<Ac_output_merge_with_model_format.tab> --plekmodel=<Ac_output_PLEK_model_Dic_tr4_ml200_k5_balanced.txt> --output=<outputfile.tab>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print("This script has 3 functions the first function take the annotation file from blastx and interpro combine them with the prediction result from CPC2 and PLEK(without the model).\nThe second function merging_table_ouput_with_PLEK_model_output combine the result of the last function after a formating of the columns (space to tab) with the prediction results from PLEK with the model. The third function take the result of second and add the result of the PLEK prediction with a balanced model.\nYou can choose the function to use with the option merging like followed :")
            print('merge_annot.py --merging="without_model" --tr_inputfile=<tr_label.fasta> --blastx=<non_annot_Ac_seqid.list> --interpro=<Ac_interpro_TSV.output.condense.list> --cpc2=<tr_label_cpc2_r.tab> --plek=<Ac_output_PLEK_tr4_ml200_k5.txt> --output=<outputfile.tab>')
            print('merge_annot.py --merging="with_model" --mergedtab=<Ac_output_merge_without_model_format.tab> --plekmodel=<Ac_output_PLEK_model_Dic_tr4_ml200_k5.txt> --output=<outputfile.tab>')
            print('merge_annot.py --merging="with_model_balanced" --mergedtab=<Ac_output_merge_with_model_format.tab> --plekmodel=<Ac_output_PLEK_model_Dic_tr4_ml200_k5_balanced.txt> --output=<outputfile.tab>')
            sys.exit()
        elif opt in ( "--merging") :
            option = arg
        elif opt in ( "--mergedtab"):
            mergedtab = arg
        elif opt in ( "--tr_inputfile"):
            tr_inputfile = arg
        elif opt in ( "--blastx"):
            blastx = arg
        elif opt in ( "--cpc2"):
            cpc2 = arg
        elif opt in ( "--plek"):
            plek = arg
        elif opt in ( "--plekmodel"):
            plekmodel = arg
        elif opt in ("-o", "--output"):
            output = arg


    if (option == "without_model") :
        merging_table(tr_inputfile, blastx, interpro, cpc2, plek, output)
    elif (option == "with_model") :
        merging_table_ouput_with_PLEK_model_output(mergedtab, plekmodel, output)
    elif (option == "with_model_balanced") :
        merging_table_ouput_with_PLEK_model_balanced_output(mergedtab, plekmodel, output)
    else :
        print("you need to chose which merging function you need you need")
        sys.exit(2)

if __name__ == "__main__":
   main(sys.argv[1:])




#merging_table("01.LABELED_DATA/Ac_labeled/Ac_tr_label.fasta","03.NON_ANNOTATED_DATA/Non_annotated_Ac/non_annot_Ac_seqid.list","04.LNCRNA_PREDICTION/4.3.results_observations/Ac/Ac_interpro_TSV.output.condense.list","04.LNCRNA_PREDICTION/4.1.cpc2_pred/Ac/Ac_tr_label_cpc2_r.tab","04.LNCRNA_PREDICTION/4.2.plek_pred/Ac/Ac_output_PLEK_tr4_ml200_k5.txt" ,"04.LNCRNA_PREDICTION/4.3.results_observations/Ac/Ac_output_merge.tab")
#merging_table("01.LABELED_DATA/AcPa_labeled/AcPa_tr_label.fasta","03.NON_ANNOTATED_DATA/Non_annotated_AcPa/non_annot_AcPa_seqid.list","04.LNCRNA_PREDICTION/4.3.results_observations/AcPa/AcPa_interpro_TSV.output.condense.list","04.LNCRNA_PREDICTION/4.1.cpc2_pred/AcPa/AcPa_tr_label_cpc2_r.tab","04.LNCRNA_PREDICTION/4.2.plek_pred/AcPa/AcPa_output_PLEK_tr4_ml200_k5.txt" ,"04.LNCRNA_PREDICTION/4.3.results_observations/AcPa/AcPa_output_merge.tab")

#merging_table("01.LABELED_DATA/Ac_labeled/Ac_tr_label.fasta","03.NON_ANNOTATED_DATA/Non_annotated_Ac/non_annot_Ac_seqid.list","04.LNCRNA_PREDICTION/4.3.results_observations/Ac/Ac_interpro_TSV.output.condense.list","04.LNCRNA_PREDICTION/4.1.cpc2_pred/Ac/Ac_tr_label_cpc2_r.tab","04.LNCRNA_PREDICTION/4.2.plek_pred/with_model/Ac/Ac_output_PLEK_model_Dic_tr4_ml200_k5.txt" ,"04.LNCRNA_PREDICTION/4.3.results_observations/Ac/Ac_output_merge_model.tab")
#merging_table("01.LABELED_DATA/AcPa_labeled/AcPa_tr_label.fasta","03.NON_ANNOTATED_DATA/Non_annotated_AcPa/non_annot_AcPa_seqid.list","04.LNCRNA_PREDICTION/4.3.results_observations/AcPa/AcPa_interpro_TSV.output.condense.list","04.LNCRNA_PREDICTION/4.1.cpc2_pred/AcPa/AcPa_tr_label_cpc2_r.tab","04.LNCRNA_PREDICTION/4.2.plek_pred/with_model/AcPa/AcPa_output_PLEK_model_Dic_tr4_ml200_k5.txt" ,"04.LNCRNA_PREDICTION/4.3.results_observations/AcPa/AcPa_output_merge_model.tab")

#merging_table_output_with_PLEK_model_output("04.LNCRNA_PREDICTION/4.3.results_observations/without_model/Ac/Ac_output_merge_without_model_format.tab", "04.LNCRNA_PREDICTION/4.2.plek_pred/with_model/Ac/Ac_output_PLEK_model_Dic_tr4_ml200_k5.txt", "04.LNCRNA_PREDICTION/4.3.results_observations/with_model/Ac/Ac_output_merge_model.tab")


#merging_table_output_with_PLEK_model_output("04.LNCRNA_PREDICTION/4.3.results_observations/without_model/AcPa/AcPa_output_merge_without_model_format.tab", "04.LNCRNA_PREDICTION/4.2.plek_pred/with_model/AcPa/AcPa_output_PLEK_model_Dic_tr4_ml200_k5.txt", "04.LNCRNA_PREDICTION/4.3.results_observations/with_model/AcPa/AcPa_output_merge_model.tab")

