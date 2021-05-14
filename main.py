#main.py

from core import app

from Bio import SeqIO
import pandas as pd
import re
import importlib

def fasta2handle(file_path):
    with open(file_path) as handle:
        return list(SeqIO.parse(handle, "fasta"))

def gff_maker2handle(file_path):
    with open(file_path) as handle:
        gff = pd.read_csv(handle, sep="\t")
        gff = gff.reset_index()
        gff.columns = ["seq_id", "source", "ft_type", "start", "stop", "score", "strand", "phase", "attr"]
        gff = gff.sort_values(by=["seq_id"]).drop(["source", "score", "phase"], axis=1).dropna()
        gff["sub_seq_id"] = [re.split(r':',re.search("^ID=.*?;", x).group(0)[3:-1])[0] for x in gff["attr"]]
        return gff[["seq_id", "sub_seq_id", "ft_type", "start", "stop", "strand"]].sort_values(by=["seq_id", "sub_seq_id", "ft_type"]).set_index(["seq_id", "sub_seq_id", "ft_type", "start", "stop", "strand"]) 
    
def tab_pannzer2handle(file_path):
    with open(file_path) as handle:
        anno = pd.read_csv(handle, sep="\t")
        anno = anno.sort_values(by=["qpid", "type"]).drop(["score", "PPV"], axis=1).set_index(["qpid", "type", "id", "desc"])
        return anno

if __name__ == "__main__":

    _GLOBALS_ = {
        "handles":dict(),
        "plugins":dict(),
        "metadata":dict()
    }

    _GLOBALS_["handles"]["fasta"] = fasta2handle("files/sequences.fasta")
    _GLOBALS_["handles"]["gff_maker"] = gff_maker2handle("files/data.gff")
    _GLOBALS_["handles"]["tab_panzer"] = tab_pannzer2handle("files/anno.out")

    _GLOBALS_["plugins"]["read_fasta"] = importlib.import_module(".read_fasta","Plugins").Plugin()

    _GLOBALS_["plugins"]["read_gff_maker_3UTR"] = importlib.import_module(".read_gff_maker_3UTR","Plugins").Plugin()
    _GLOBALS_["plugins"]["read_gff_maker_5UTR"] = importlib.import_module(".read_gff_maker_5UTR","Plugins").Plugin()
    _GLOBALS_["plugins"]["read_gff_maker_CDS"] = importlib.import_module(".read_gff_maker_CDS","Plugins").Plugin()
    _GLOBALS_["plugins"]["read_gff_maker_exon"] = importlib.import_module(".read_gff_maker_exon","Plugins").Plugin()
    _GLOBALS_["plugins"]["read_gff_maker_gene"] = importlib.import_module(".read_gff_maker_gene","Plugins").Plugin()
    _GLOBALS_["plugins"]["read_gff_maker_main"] = importlib.import_module(".read_gff_maker_main","Plugins").Plugin()
    _GLOBALS_["plugins"]["read_gff_maker_mRNA"] = importlib.import_module(".read_gff_maker_mRNA","Plugins").Plugin()

    _GLOBALS_["plugins"]["read_tab_pannzer_CDS"] = importlib.import_module(".read_tab_pannzer_CDS","Plugins").Plugin()
    _GLOBALS_["plugins"]["read_tab_pannzer_gene"] = importlib.import_module(".read_tab_pannzer_gene","Plugins").Plugin()


    _GLOBALS_["metadata"]["project"] = "temp"
    _GLOBALS_["metadata"]["division"] = "INV"
    _GLOBALS_["metadata"]["taxonomy"] = "29031"
    _GLOBALS_["metadata"]["organism"] = "Phlebotomus papatasi"
    _GLOBALS_["metadata"]["molecule_type"] = "genomic DNA"
    _GLOBALS_["metadata"]["topology"] = "linear"
    _GLOBALS_["metadata"]["description"] = "description"
    _GLOBALS_["metadata"]["transl_table"] = 0

    app = app(
        [
            (_GLOBALS_["plugins"]["read_fasta"], _GLOBALS_["handles"]["fasta"], _GLOBALS_["metadata"], [
                (_GLOBALS_["plugins"]["read_gff_maker_main"], _GLOBALS_["handles"]["gff_maker"], _GLOBALS_["metadata"], [
                    (_GLOBALS_["plugins"]["read_gff_maker_gene"], _GLOBALS_["handles"]["gff_maker"], _GLOBALS_["metadata"], [
                        (_GLOBALS_["plugins"]["read_tab_pannzer_gene"], _GLOBALS_["handles"]["tab_panzer"], _GLOBALS_["metadata"], [])
                    ]),
                    (_GLOBALS_["plugins"]["read_gff_maker_mRNA"], _GLOBALS_["handles"]["gff_maker"], _GLOBALS_["metadata"], []),
                    (_GLOBALS_["plugins"]["read_gff_maker_CDS"], _GLOBALS_["handles"]["gff_maker"], _GLOBALS_["metadata"], [
                        (_GLOBALS_["plugins"]["read_tab_pannzer_CDS"], _GLOBALS_["handles"]["tab_panzer"], _GLOBALS_["metadata"], [])
                    ]),
                    (_GLOBALS_["plugins"]["read_gff_maker_3UTR"], _GLOBALS_["handles"]["gff_maker"], _GLOBALS_["metadata"], []),
                    (_GLOBALS_["plugins"]["read_gff_maker_5UTR"], _GLOBALS_["handles"]["gff_maker"], _GLOBALS_["metadata"], [])
                ]),
                (_GLOBALS_["plugins"]["read_gff_maker_exon"], _GLOBALS_["handles"]["gff_maker"], _GLOBALS_["metadata"], [])
            ])
        ])
    app.run()