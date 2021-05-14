#fasta2handle.py

from Bio import SeqIO

class Plugin:

    def process(self, file_path):
        with open(file_path) as handle:
            return list(SeqIO.parse(handle, "fasta"))