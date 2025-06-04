from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO, AlignIO
from Bio.Align import PairwiseAligner
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

human_cyt_c = SeqIO.read("P99999.fasta", format="fasta") #import protein
result_handle = NCBIWWW.qblast("blastp", "nr", human_cyt_c.format("fasta"), hitlist_size=8) #blast gia homologous
blast_records = NCBIXML.read(result_handle)
homologs = []

for alignment in blast_records.alignments:
    for hsp in alignment.hsps:
        seq_record = SeqRecord(Seq(hsp.sbjct),
                               id=alignment.hit_id,
                               description=alignment.hit_def)
        homologs.append(seq_record)
        break
    if len(homologs) >= 8: #8 omologa
        break


homologs.insert(0, human_cyt_c)
SeqIO.write(homologs, "cytc_homologs.fasta", "fasta")# akolou8ies se fasta

clustalo_cline = ClustalOmegaCommandline(cmd=r"C:\Users\Paris\Downloads\clustal-omega-1.2.2-win64\clustal-omega-1.2.2-win64\clustalo.exe", #crustalo
                                         infile="cytc_homologs.fasta",
                                         outfile="cytc_homologs.aln5",
                                         verbose=True,
                                         auto=True,
                                         outfmt="clu")  # clustal format
stdout, stderr = clustalo_cline()


alignment = AlignIO.read("cytc_homologs.aln5", "clustal") #load allignment

aligner = PairwiseAligner() #pairs
aligner.mode = 'global'

for i in range(len(homologs)): #write em files
    for j in range(i + 1, len(homologs)):
        aln = aligner.align(homologs[i].seq, homologs[j].seq)[0]
        with open(f"pairwise_{i}_{j}.aln5", "w") as f:
            f.write(f"> {homologs[i].id}\n{aln.aligned[0]}\n")
            f.write(f"> {homologs[j].id}\n{aln.aligned[1]}\n")
            f.write(f"\nScore: {aln.score}\n")
