import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import glob, os, sys, getopt

def merging_table(sRNA_tab, filecpc2, fileplek, fileplekmodel_non_balanced, plekmodel_balanced , output):
    #opening the output
    out=open(output,"w")
    
    
    #creation of the dataframe with only the id of the seq and its length
    df_input = pd.read_table(sRNA_tab,dtype={'Seq_Length': str},names=['line','ID','Seq_Length','Annotation_Blastx', 'Annotation_Interpro_TSV','Sequence'])
    df_input = df_input.drop(columns="line")
    df_input = df_input.iloc[1:]
    print("Analysis of the annotation files starting ...")

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
    print("Plek without model analysis starting.\nCreation of the dataframe commencing ...")

    plek = open(fileplek, 'r')
    pleklist=list()
    for lineplek in plek.readlines():
        fullidplek =  lineplek.split('\t')[2]
        idfastaplek =  fullidplek.split(' ')[0]
        idplek = idfastaplek.split('>')[1]
        pleklist.append([idplek,lineplek.split('\t')[1],lineplek.split('\t')[0]])
    df_plek= pd.DataFrame(pleklist,columns=['ID','PLEK_score','PLEK_prediction'])
    


    plek_model_non_balanced = open(fileplekmodel_non_balanced, 'r')
    plek_model_non_balanced_list=list()
    for lineplek in plek_model_non_balanced.readlines():
        fullidplek =  lineplek.split('\t')[2]
        idfastaplek =  fullidplek.split(' ')[0]
        idplek = idfastaplek.split('>')[1]
        plek_model_non_balanced_list.append([idplek,lineplek.split('\t')[1],lineplek.split('\t')[0]])
    df_plek_model_non_balanced= pd.DataFrame(plek_model_non_balanced_list,columns=['ID','PLEK_score_model','PLEK_prediction_model'])
    

    
    plek_model_balanced = open(plekmodel_balanced, 'r')
    plek_model_balancedlist=list()
    for lineplek in plek_model_balanced.readlines():
        fullidplek =  lineplek.split('\t')[2]
        idfastaplek =  fullidplek.split(' ')[0]
        idplek = idfastaplek.split('>')[1]
        plek_model_balancedlist.append([idplek,lineplek.split('\t')[1],lineplek.split('\t')[0]])
    df_plek_model_balanced= pd.DataFrame(plek_model_balancedlist,columns=['ID','PLEK_score_model_balanced','PLEK_prediction_model_balanced'])
    

    #merging all the dataframe into one
    print("Dataframe merging started ...")
    df_output = pd.merge(df_input,df_cpc, how='outer')
    df_output = pd.merge(df_output,df_plek, how='outer')
    df_output = pd.merge(df_output,df_plek_model_non_balanced, how='outer')
    df_output = pd.merge(df_output,df_plek_model_balanced, how='outer')
    print("Merging finished ...")
    
    print("Preview of the merged dataframe ...")
    print(df_output)
    
    #writing of the output file 
    out.write(df_output.to_string())


#calling the function and its parameters
def main(argv):

    sRNA_tab=''
    cpc2=''
    plek=''
    plekmodel_non_balanced=''
    plekmodel_balanced=''
    output=''
	
    try:
        opts, args = getopt.getopt(argv,"o",["sRNA_inputfile=","cpc2=", "plek=", "plekmodel_non_balanced=", "plekmodel_balanced=", "output="])
		
    except getopt.GetoptError: 
        print('merge_annot.py --sRNA_inputfile --cpc2 --plek --plekmodel_non_balanced --plekmodel_balanced --output')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('merge_annot.py --sRNA_inputfile --cpc2 --plek --plekmodel_non_balanced --plekmodel_balanced --output')
            sys.exit()
        
        elif opt in ( "--sRNA_inputfile"):
            sRNA_tab = arg
        elif opt in ( "--cpc2"):
            cpc2 = arg
        elif opt in ( "--plek"):
            plek = arg
        elif opt in ( "--plekmodel_non_balanced"):
            plekmodel_non_balanced = arg
        elif opt in ( "--plekmodel_balanced"):
            plekmodel_balanced = arg
        elif opt in ("-o", "--output"):
            output = arg


    merging_table(sRNA_tab, cpc2, plek, plekmodel_non_balanced, plekmodel_balanced , output)
    
if __name__ == "__main__":
   main(sys.argv[1:])
