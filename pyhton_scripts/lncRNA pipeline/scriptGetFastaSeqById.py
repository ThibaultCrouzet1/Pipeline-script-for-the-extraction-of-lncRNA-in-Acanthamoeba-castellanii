import sys, getopt
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import glob, os


def filter_seq(fastaseq, idseqfile, output_file):
    #must give the path to the fasta file
    #"Ac/Transcriptome_Ac.fasta"
    print("Openning the fasta file ...")
    records = list(SeqIO.parse( fastaseq, "fasta"))

    print("openning the id file ...")
    idfile = open(idseqfile,"r")

    
    #creation of the output file that will containt all the annoted sequences
    output=open(output_file,"w")


    #search and indentification of the sequences by their id
    for line in idfile.readlines() :
        seqid = line.split(' ')[0].split('\n')[0]
        for record in records:
            if record.id == seqid:
                print(line)
                print(record)
                SeqIO.write(record, output, "fasta")
                break
    
    
    #closing the opened file
    idfile.close()
    output.close() 

#calling the function and its parameters
def main(argv):
    tr_inputfile = ''
    idlist_inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"t:i:o:",["tr_inputfile=","idlist_inputfile=","outputfile="])
    except getopt.GetoptError:
        print('scriptGetFastaSeqById.py -t <tr_inputfile.fasta> -i <idlist_inputfile> -o <outputfile>')
        print('scriptGetFastaSeqById.py --tr_inputfile=<tr_inputfile.fasta> --idlist_inputfile=<idlist_inputfile> --outputfile=<outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('scriptGetFastaSeqById.py -t <tr_inputfile.fasta> -i <idlist_inputfile> -o <outputfile>')
            print('scriptGetFastaSeqById.py --tr_inputfile <tr_inputfile.fasta> --idlist_inputfile <idlist_inputfile> --outputfile <outputfile>')
            sys.exit()
        elif opt in ("-t", "--tr_inputfile"):
            tr_inputfile = arg
        elif opt in ("-i", "--idlist_inputfile"):
            idlist_inputfile = arg
        elif opt in ("-o", "--outputfile"):
            outputfile = arg

    filter_seq(tr_inputfile,idlist_inputfile,outputfile)

if __name__ == "__main__":
   main(sys.argv[1:])

#example of parameters
#"../00.RAW_DATA/Ac/Ac_tr.fasta", "../03.NON_ANNOTATED_DATA/NON_ANNOTATED_AC_ID_SEQ.list", "../03.NON_ANNOTATED_DATA/fastaseqAC.fasta"
#"../00.RAW_DATA/AcPa/AcPa_tr.fasta", "../03.NON_ANNOTATED_DATA/NON_ANNOTATED_ACPA_ID_SEQ.list", "../03.NON_ANNOTATED_DATA/fastaseqACPA.fasta"