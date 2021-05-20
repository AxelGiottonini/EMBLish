#fasta2handle.py

from Bio import SeqIO
from Plugins.__plugin__ import __Plugin__

class Plugin(__Plugin__):

    def process(self, file_path):
        with open(file_path) as handle:
            return list(SeqIO.parse(handle, "fasta"))