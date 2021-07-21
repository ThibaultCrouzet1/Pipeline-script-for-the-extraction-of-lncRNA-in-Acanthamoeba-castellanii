import sys, getopt
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import glob, os

def seqInBothTr(tr_Ac, tr_AcPa, output_file):
    #must give the path to the fasta file
    #"Ac/Transcriptome_Ac.fasta"
    print("Openning the fasta file ...")
    Ac_records = list(SeqIO.parse( tr_Ac, "fasta"))
    AcPa_records = list(SeqIO.parse( tr_AcPa, "fasta"))

        
    #creation of the output file that will containt all the annoted sequences
    output=open(output_file,"w")


    #search and indentification of the sequences by their id
    print("extraction of the seuqence that are in both expressed transcript ...")
    for Ac_r in Ac_records :
        for AcPa_r in AcPa_records:
            if Ac_r.seq == AcPa_r.seq:
                print(Ac_r.id)
                SeqIO.write(Ac_r, output, "fasta")
                
    
    
    #closing the opened file
    output.close() 

#calling the function and its parameters
def main(argv):
    tr_Ac = ''
    tr_AcPa = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"a:p:o:",["tr_Ac=","tr_AcPa=","outputfile="])
    except getopt.GetoptError:
        print('scriptGetFastaSeqById.py -a <tr_Ac.fasta> -p <tr_AcPa.fasta> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('scriptGetFastaSeqById.py -a <tr_Ac.fasta> -p <tr_AcPa.fasta> -o <outputfile>')
            sys.exit()
        elif opt in ("-a", "--tr_Ac"):
            tr_Ac = arg
        elif opt in ("-p", "--tr_AcPa"):
            tr_AcPa = arg
        elif opt in ("-o", "--outputfile"):
            outputfile = arg

    seqInBothTr(tr_Ac,tr_AcPa,outputfile)

if __name__ == "__main__":
   main(sys.argv[1:])