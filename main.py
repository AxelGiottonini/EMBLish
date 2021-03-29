#!/usr/bin/python3
# -*-coding:utf8

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqRecord import SeqRecord 
from joblib import Parallel, delayed
from tqdm import tqdm
import multiprocessing
import os
import pandas as pd
import re
import sys

def read_anno(file):
    anno = pd.read_csv(file, sep="\t")
    anno = anno.sort_values(by=["qpid", "type"]).drop(["score", "PPV"], axis=1).set_index(["qpid", "type", "id", "desc"])
    return anno

def read_gff(file):
    gff = pd.read_csv(file, sep="\t")
    gff = gff.reset_index()
    gff.columns = ["seq_id", "source", "ft_type", "start", "stop", "score", "strand", "phase", "attr"]
    gff = gff.sort_values(by=["seq_id"]).drop(["source", "score", "phase"], axis=1).dropna()
    gff["sub_seq_id"] = [re.split(r':',re.search("^ID=.*?;", x).group(0)[3:-1])[0] for x in gff["attr"]]
    gff = gff[["seq_id", "sub_seq_id", "ft_type", "start", "stop", "strand"]].sort_values(by=["seq_id", "sub_seq_id", "ft_type"]).set_index(["seq_id", "sub_seq_id", "ft_type", "start", "stop", "strand"]) 
    return gff

def read_config(file):
    config = {"description":None, "division":None, "molecule_type":None, "organism":None,  "project":None, "taxonomy":None, "topology":None, "transl_table":None}
    for line in [line.rsplit("\n")[0] for line in file.readlines()]:
        config[re.split(r":", line)[0].lower()]=re.split(r":", line)[1]
    return config

def mergeLocations(_locationArray_):
    return _locationArray_[0] if len(_locationArray_) == 1 else CompoundLocation(_locationArray_)

def init_features(sec, gene, gff):
    ft_table = gff.loc[(sec, gene),:].reset_index()
    ft_table_mRNA = ft_table[ft_table["ft_type"] == "mRNA"].iloc[0,:]
    ft_table_CDSs = [x for i, x in ft_table[ft_table["ft_type"] == "CDS"].iterrows()]
    tmp_ft_table_3UTR = ft_table[ft_table["ft_type"] == "three_prime_UTR"]
    ft_table_3UTR = pd.Series(dtype="float64") if tmp_ft_table_3UTR.empty else tmp_ft_table_3UTR.iloc[0,:]
    tmp_ft_table_5UTR = ft_table[ft_table["ft_type"] == "five_prime_UTR"]
    ft_table_5UTR = pd.Series(dtype="float64") if tmp_ft_table_5UTR.empty else tmp_ft_table_5UTR.iloc[0,:]
    
    return {
        "location":FeatureLocation(int(ft_table_mRNA[1]),int(ft_table_mRNA[2]),(1,-1)[ft_table_mRNA[3] == "-"]), 
        "qualifiers":{"gene":gene,"note":list()},
        "type":"gene"
    },{
        "location":mergeLocations([FeatureLocation(int(CDS[1]),int(CDS[2]),(1,-1)[CDS[3] == "-"]) for CDS in ft_table_CDSs]), 
        "qualifiers":{"gene":gene}, 
        "type":"mRNA"
    },{
        "location":mergeLocations([FeatureLocation(int(CDS[1]),int(CDS[2]),(1,-1)[CDS[3] == "-"]) for CDS in ft_table_CDSs]), 
        "qualifiers":{"gene":gene,"product":list(),"note":list(),"db_xref":list(),"translation":list(),"transl_table":11}, 
        "type":"CDS"
    },{
        "location": None if ft_table_3UTR.empty else FeatureLocation(int(ft_table_3UTR[1]), int(ft_table_3UTR[2]), (1,-1)[ft_table_3UTR[3] == "-"]),
        "qualifiers":{"gene":gene}, 
        "type":"3'UTR"
    },{
        "location":None if ft_table_5UTR.empty else FeatureLocation(int(ft_table_5UTR[1]), int(ft_table_5UTR[2]), (1,-1)[ft_table_5UTR[3] == "-"]),
        "qualifiers":{"gene":gene}, 
        "type":"5'UTR"}

def merge(record, anno, gff, conf, out_dir):
    sec = record.id
    
    #Record initialisation
    _record_ = SeqRecord(
        record.seq,
        id=sec,
        dbxrefs=["Project:" + conf["project"]],
        annotations={"division":conf["division"],"molecule_type":conf["molecule_type"],"organism":conf["organism"],"taxonomy":conf["taxonomy"],"topology":conf["topology"]},
        description=conf["description"])
    
    #Source feature
    ft_table = gff.loc[(sec, slice(None), "contig"), :].reset_index()
    _source_ = SeqFeature(FeatureLocation(int(ft_table.iloc[0,3]),int(ft_table.iloc[0,4]),(1,-1)[ft_table.iloc[0,5] == "-"]),type="source",qualifiers={"organism":"test","mol_type":"genomic DNA","db_xref":list()})
    _record_.features.append(_source_)
    
    #GENE/MRNA/CDS/3UTR/5UTR features
    for gene in gff.loc[(sec, slice(None), "gene"),:].reset_index()["sub_seq_id"].apply(lambda x: x+"-mRNA-1"):
        _gene_, _mRNA_, _CDS_, _3UTR_, _5UTR_ = init_features(sec, gene, gff)
        
        try:
            anno_table = anno.loc[(gene),:].reset_index()
            anno_bp = anno.loc[(gene, "BP_ARGOT"),:].reset_index()["id"]
            anno_cc = anno.loc[(gene, "CC_ARGOT"),:].reset_index()["id"]
            anno_mf = anno.loc[(gene, "MF_ARGOT"),:].reset_index()["id"]
            _CDS_["qualifiers"]["db_xref"] = ["GO:" + str(go) for go in pd.concat([anno_bp, anno_cc, anno_mf])]
            _CDS_["qualifiers"]["translation"] = anno.loc[(gene, "qseq"),:].reset_index().iloc[0,1]
            _CDS_["qualifiers"]["transl_table"] = conf["transl_table"]
            _gene_["qualifiers"]["note"] = _CDS_["qualifiers"]["product"] = anno.loc[(gene, "DE"),:].reset_index().iloc[0,1]
        except KeyError:
            pass
        
        for feature in [_gene_, _mRNA_, _CDS_, _3UTR_, _5UTR_]:
            if feature["location"]: 
                _record_.features.append(SeqFeature(feature["location"], type=feature["type"], qualifiers=feature["qualifiers"]))
    
    #Exon feature
    ft_table = gff.loc[sec, slice(None), "exon"].reset_index()
    _record_.features.extend([SeqFeature(FeatureLocation(int(exon["start"]), int(exon["stop"]), (1,-1)[exon["strand"]=="-"]), type="exon", qualifiers={}) for i,exon in ft_table.iterrows()])
     
    #Print EMBL entry in output folder
    with open(out_dir + "/" + sec + ".dat", "w") as file:
        print(_record_.format("embl"), file=file)
        file.close()   
    
if __name__ == "__main__":

    GFF_FILE = FASTA_FILE = ANNO_FILE = CONF_FILE = None
    OUT_DIR = "out"

    args = sys.argv[1:]
    for i in [0,2,4,6,8]:
        if args[i] in ["-gff", "-g"]: GFF_FILE = args[i+1]
        elif args[i] in ["-fasta", "-f"]: FASTA_FILE = args[i+1]
        elif args[i] in ["-anno", "-a"]: ANNO_FILE = args[i+1]
        elif args[i] in ["-conf", "-c"]: CONF_FILE = args[i+1]
        elif args[i] in ["-out", "-o"]: OUT_DIR = args[i+1]

    if not os.path.exists(OUT_DIR):
        os.makedirs(OUT_DIR)

    CONF = None
    with open(CONF_FILE) as conf_file:
        CONF = read_config(conf_file)
        conf_file.close()
    print("Configuration file reading : DONE!")

    ANNO = None
    with open(ANNO_FILE) as anno_file:
        ANNO = read_anno(anno_file)
        anno_file.close()
    print("Annotation file reading : DONE!")

    GFF = None
    with open(GFF_FILE) as gff_file:
        GFF = read_gff(gff_file)
        gff_file.close()
    print("Prediction file reading : DONE!")

    with open(FASTA_FILE) as fasta_file:
        num_cores = multiprocessing.cpu_count()
        records = list(SeqIO.parse(fasta_file, "fasta"))
        processed_list = Parallel(n_jobs=num_cores)(delayed(merge)(record, ANNO, GFF, CONF, OUT_DIR) for record in tqdm(records))
        fasta_file.close()