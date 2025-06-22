from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import re

mystery_seq = """>mystery_sequence
GATGAPGIAGAPGFPGARGAPGPQGPSGAPGPKXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXGVQGPPGPQGPR
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXGSAGPPGATGFP
GAAGRXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXGVVGLPGQR"""


sequence = "".join(mystery_seq.split('\n')[1:]).replace(' ', '').replace('\r', '')
result_handle = NCBIWWW.qblast("blastp", "nr", sequence)

with open("blast_result.xml", "w") as out_handle:
    out_handle.write(result_handle.read())
result_handle.close()

with open("blast_result.xml") as result_handle:
    blast_records = NCBIXML.parse(result_handle)
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            match = re.search(r'\[(.*?)\]', alignment.hit_def)
            if match:
                organism_name = match.group(1)
                if organism_name.startswith("Tyr"):
                    for hsp in alignment.hsps:
                        print(f"Organism: {organism_name}")
                        print(f"Protein: {alignment.hit_def}")
                        print(f"Score: {hsp.score}")
                        print(f"E-value: {hsp.expect}")
                        print(f"Subject Sequence (part): {hsp.sbjct[:60]}...")
