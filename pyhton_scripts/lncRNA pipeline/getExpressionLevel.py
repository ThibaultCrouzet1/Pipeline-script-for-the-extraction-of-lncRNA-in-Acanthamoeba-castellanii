from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import glob, os, sys, getopt

#the function that will extract the expression level of the lncRNA sequences
def getExpressionLevel(fasta,expressionlevel,output):
    #opening the fasta file
    records = list(SeqIO.parse(fasta, "fasta"))
    
    #opening the expression level file
    expre = open(expressionlevel,'r')

    #opening the output
    out = open(output,"w")

    #getting the first line to have the description of the columns
    out.write(expre.readline())

    #filtering
    print("Extracting ...")
    for line in expre.readlines() :
        for record in records :
            if record.id in line :
                out.write(line)
                break
    print("Done")
    expre.close()
    out.close()




#calling the function and its parameters
def main(argv):
    
    fasta=''
    expreLevel=''
    output=''
	
    try:
        opts, args = getopt.getopt(argv,'f:e:o:',["fasta=","expressionlevel=","output="])
		
    except getopt.GetoptError: 
        print('getExpressionLevel.py -f <lncRNA.fasta> -e <expressionLevel_inputfile.txt/tab> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('getExpressionLevel.py -f <lncRNA.fasta> -e <expressionLevel_inputfile.txt/tab> -o <outputfile>')
            sys.exit()
        elif opt in ( "-f","--fasta") :
            fasta = arg
        elif opt in ( "-e","--expressionlevel"):
            expreLevel = arg
        elif opt in ("-o","--output"):
            output = arg

    getExpressionLevel(fasta,expreLevel,output)


if __name__ == "__main__":
   main(sys.argv[1:])


#getExpressionLevel.py -f /home/crouzet/05.LNCRNA_EXTRACTION/5.2.KNOWN_LNCRNA/Ac/Ac_lncRNA_extracted_balanced_seq.fa -e /home/crouzet/06.EXPRESSION_LEVEL/6.2.TRANSCRIPTS_LABELED_EXPRESSION_LEVEL/Ac/Ac_labeled_diffExpr.P1e-5_C2.matrix.log2.centered.txt -o /home/crouzet/06.EXPRESSION_LEVEL/6.3.LNCRNA_EXPRESSION_LEVEL/Ac/Ac_lncRNA_expression_level.tab
#getExpressionLevel.py -f /home/crouzet/05.LNCRNA_EXTRACTION/5.2.KNOWN_LNCRNA/AcPa/AcPa_lncRNA_extracted_balanced_seq.fa -e /home/crouzet/06.EXPRESSION_LEVEL/6.2.TRANSCRIPTS_LABELED_EXPRESSION_LEVEL/AcPa/AcPa_labeled_diffExpr.P1e-5_C2.matrix.log2.centered.txt -o /home/crouzet/06.EXPRESSION_LEVEL/6.3.LNCRNA_EXPRESSION_LEVEL/AcPa/AcPa_lncRNA_expression_level.tab
