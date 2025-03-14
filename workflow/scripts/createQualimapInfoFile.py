import sys
# import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import Bio.SeqUtils as SeqUtils
import HTSeq
import numpy as np


def idsContainGiven(givenId, transcriptIds):
    for tId in transcriptIds:
        if givenId.find(tId) != -1:
            return True
    return False

if __name__ == "__main__":
    
    descriptionText = "The script extracts features from a GTF file and a FASTA file into Qualimap annotation format. Note: exons have to be sorted according to exon number! This important for correct reverse transcribed cDNA sequences extraction."

    # parser = argparse.ArgumentParser(description = descriptionText,formatter_class=argparse.RawDescriptionHelpFormatter)
    
    # parser.add_argument("-g", action="store", required="true", dest="gtfFile",
    #     help="Input file with list of genes in GTF format")
    # parser.add_argument("-f", action="store", required="true", dest="fastaFile",
    #     help="Input genome sequence. ")
    # parser.add_argument("-o", action="store", dest="outFile", default="annotations.txt",
    #     help="Output file. Default is annotations.txt")

    # parser.add_argument("--filter", action="store", dest="filterStr", default="",
    #                     help="Comma-separted list of entries to filter from GTF file \
    #                     based on given attribute id")

    # parser.add_argument("--ignore-strange-chrom", action="store_true", default=False,
    #     dest="ignoreStrangeChromosomes", help="All chromosomes except numbered and X,Y,MT are ignored ")
    
    # args = parser.parse_args()
    
    # print(args)
    sys.stdout = open(snakemake.log[0], 'w')
    sys.stderr = sys.stdout

    gtfFileName = snakemake.input[0]
    fastaFileName = snakemake.input[1]
    outFileName = snakemake.output[0]
    attr_id = "gene_id"

    
    # parse GTF file

    gtf_file = HTSeq.GFF_Reader( gtfFileName )

    features = {}

    filtered_transcripts = ""
    if filtered_transcripts:
        print(f"Filtering for: {filtered_transcripts}")

    for feature in gtf_file:
        if feature.type == 'exon':
            geneName = feature.attr[ attr_id ]
            #print transcriptName
            if geneName in features:
                features[geneName].append(feature)
            else:
                features[geneName] = [feature]

    # load & save sequences 

    seqData = SeqIO.to_dict(SeqIO.parse(fastaFileName, "fasta"))
    
    outFile = open(outFileName, "w")

    header = "\"%s\"\t\"%s\"\t\"%s\"\n" % ("biotypes","length","gc")
    outFile.write(header)


    for geneId in features:
       
        exons = features[geneId]

        print(f"Processing {geneId}")

        if len(exons) == 0:
            continue

        biotype = "unknown"
        length = 0
        transcripts = {}

        for exon in exons:
            transcriptId = exon.attr["transcript_id"]
            
            tSeq = transcripts.get(transcriptId, Seq(""))

            iv = exon.iv
            seqName = iv.chrom
            if seqName in seqData:
                #print "Exon (%s,%d,%d) " % (iv.chrom,iv.start,iv.end)
                buf = seqData[ iv.chrom ].seq[ iv.start  : iv.end ]
                if iv.strand == '-':
                    buf = buf.reverse_complement()
                tSeq += buf
            else:
                print(f"Can not locate sequence {seqName} in {fastaFileName}, skipping region...")
            transcripts[transcriptId] = tSeq

        gc_array = []
        lengths = []

        for tSeq in transcripts.values():
            if len(tSeq) > 5 and set(tSeq.upper()) != {'N'}:
                print(tSeq)
                lengths.append( len(tSeq) )
                GC_content = SeqUtils.GC123(tSeq)[0]
                gc_array.append ( GC_content )
            

        #gene_length = len(seq_rec) / len(transcripts)
        #gc = SeqUtils.GC(seq_rec)
        
        #print gc_array, lengths
        gene_length = np.mean(lengths)
        gene_gc = np.mean(gc_array)
        
        #print len(seq_rec), len(transcripts),gc

        line = f'"{geneId}"\t"{biotype}"\t{gene_length}\t{gene_gc}\n'

        outFile.write ( line )
                   
        #outFile.flush()
        #sys.stdin.readline()

    outFile.close()


