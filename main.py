#!/usr/bin/python3
# -*-coding:utf-8

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

import os.path
import re
import subprocess
import sys


class GFF_ENTRY():
    def __init__(self, entry):
        parsed = re.split(r"\t", entry)
        self.seqid = parsed[0]
        self.type = parsed[2]
        self.start = int(parsed[3]) - 1
        self.end = int(parsed[4])
        self.strand = -1 if parsed[6] == "-" else 1
        self.attributes = parsed[8]

    def exonId(self):
        assert (self.type == "exon"), "Exon entry required"
        return ("exon_id=" + self.seqid + "-E:" + re.split(r":", re.split(r";", self.attributes)[0])[2])

    def length(self):
        return self.end - self.start

class ANNOTATION_ENTRY():
    def __init__(self, entry):
        parsed = re.split(r"\t", entry)
        self.qpid = parsed[0]
        self.type = parsed[1]
        self.id = parsed[4]
        self.desc = parsed[5]

    def keggId(self):
        assert (self.type == "EC_ARGOT"), "EC entry required"
        for id in re.split(r" ", self.id):
            if id and re.search("^EC:\d+\.\d+\.\d+\.\d+$", id): (yield "KEGG_Enzyme:" + id[3:])

def cmd(_cmd_):
    subprocess.call(_cmd_, shell=True)

def mergeLocations(_locationArray_):
    assert (len(_locationArray_) != 0), "Empty location array" 
    return _locationArray_[0] if len(_locationArray_) == 1 else CompoundLocation(_locationArray_)

def main(fastaFile, annoFile, gffFile, projectFile):
    """Script files & folder assignment and initialisation. The files have to be specified in 
    the function call, and the already existing temp & out folder will be removed with their 
    whole content."""
    SEQUENCES_FILE = fastaFile
    ANNOTATIONS_FILE = annoFile
    GFF_FILE = gffFile 
    PROJECT_FILE = projectFile   
    TEMP_FOLDER = "temp/"
    OUT_FOLDER = "out/"
    if os.path.isdir(TEMP_FOLDER): cmd("rm -r " + TEMP_FOLDER)
    if os.path.isdir(OUT_FOLDER): cmd("rm -r " + OUT_FOLDER)
    cmd("mkdir " + TEMP_FOLDER)
    cmd("mkdir " + OUT_FOLDER)
    assert os.path.isfile(SEQUENCES_FILE), "FASTA file not found!"
    assert os.path.isfile(ANNOTATIONS_FILE), "Annotations file not found!"
    assert os.path.isfile(GFF_FILE), "GFF file not found!"
    assert os.path.isfile(projectFile), "Project file not found!"
    assert os.path.isdir(TEMP_FOLDER), "Could not create temp folder!"
    assert os.path.isdir(OUT_FOLDER), "Could not create out folder!"

    """Reading of the project file and definition of the shared variables"""
    __PROJECT__ = None
    #__DATA_CLASS__ = None
    __DIVISION__ = None
    __ORGANISM__ = None
    __TAXONOMY__ = None
    __MOLECULE_TYPE__ = None
    __TOPOLOGY__ = None
    __DESCRIPTION__ = None
    with open(PROJECT_FILE, "r") as project:
        project_lines = (line.rstrip("\n") for line in project.readlines())
        try:
            while True:
                parsed = re.split(r":", next(project_lines))
                if parsed[0] == "PROJECT": __PROJECT__ = parsed[1]
                elif parsed[0] == "DIVISION": __DIVISION__ = parsed[1]
                elif parsed[0] == "TAXONOMY": __TAXONOMY__ = re.split(r"-",parsed[1])
                elif parsed[0] == "ORGANISM": __ORGANISM__ = parsed[1]
                elif parsed[0] == "MOLECULE_TYPE": __MOLECULE_TYPE__ = parsed[1]
                elif parsed[0] == "TOPOLOGY": __TOPOLOGY__ = parsed[1]
                elif parsed[0] == "DESCRIPTION": __DESCRIPTION__ = parsed[1]
        except StopIteration:
            pass
        finally:
            project.close()
    assert __PROJECT__, "Project undefined"
    assert __DIVISION__, "Division undefined"
    assert __ORGANISM__, "Organism undefined"
    assert __TAXONOMY__ , "Taxonomy undefined"
    assert __MOLECULE_TYPE__, "Molecule type undefined"
    assert __TOPOLOGY__, "Topology undefined"
    assert __DESCRIPTION__, "Description undefined"

    """Creation of the contig file containing a list of all the sequences information, 
    complementary to the fasta file"""
    CONTIG_FILE = TEMP_FOLDER + "contig.gff"
    cmd("grep \"contig\" " + GFF_FILE + " > " + CONTIG_FILE)
    assert os.path.isfile(CONTIG_FILE)

    """Line by line reading of the fasta file containing the sequence informations of the
    contigs. We should note that the sequence id has to be found in the annotations and gff
    files."""
    with open(SEQUENCES_FILE, "r") as sequences:
        sequences_lines = (line.rstrip("\n") for line in sequences.readlines())
        
        while True:
            try:
                sequence_id = next(sequences_lines)[1:]
                sequence = next(sequences_lines)

                """Creation of the files and folder required for the sequence merging and
                initialization of the record to which we already append the source feature"""
                SEQUENCE_FOLDER = TEMP_FOLDER + sequence_id + "/"
                SEQUENCE_ANNOTATIONS_FILE = SEQUENCE_FOLDER + "anno.out"
                SUB_SEQUENCES_FILE = SEQUENCE_FOLDER + "subseq.out"
                SEQUENCE_EXONS_FILE = SEQUENCE_FOLDER + "exons.gff"
                OUTPUT_FILE = OUT_FOLDER + sequence_id + ".dat"
                cmd("mkdir " + SEQUENCE_FOLDER)
                cmd("grep \"" + sequence_id + "\" " + ANNOTATIONS_FILE + " > " + SEQUENCE_ANNOTATIONS_FILE)
                cmd("grep \"original_DE\" " + SEQUENCE_ANNOTATIONS_FILE + " | cut -f1 > " + SUB_SEQUENCES_FILE)
                cmd("grep \"" + sequence_id + "\" " + GFF_FILE + " | grep \":exon:\" > " + SEQUENCE_EXONS_FILE)
                cmd("touch " + OUTPUT_FILE)
                assert os.path.isdir(SEQUENCE_FOLDER), "Could not create sequence folder!"
                assert os.path.isfile(SEQUENCE_ANNOTATIONS_FILE), "Could not create sequence annotation file!"
                assert os.path.isfile(SUB_SEQUENCES_FILE), "Could not create sub sequences file!"
                assert os.path.isfile(SEQUENCE_EXONS_FILE)
                assert os.path.isfile(OUTPUT_FILE)

                _record_ = SeqRecord(Seq(sequence), id=sequence_id)
                _record_.dbxrefs = ["Project:"+__PROJECT__]
                _record_.annotations = {
                    "molecule_type":__MOLECULE_TYPE__, 
                    "topology":__TOPOLOGY__, 
                    "data_file_division":__DIVISION__, 
                    "taxonomy":__TAXONOMY__, 
                    "organism":__ORGANISM__}
                _record_.description = __DESCRIPTION__

                source_entry = GFF_ENTRY(subprocess.check_output("grep \"" + sequence_id + "\" " + CONTIG_FILE, shell=True).decode("utf-8"))
                _source_ = SeqFeature(FeatureLocation(source_entry.start, source_entry.end, strand=source_entry.strand), type="source")
                _source_QUALIFIERS = {"organism":[__ORGANISM__], "mol_type":[__MOLECULE_TYPE__], "db_xref":list()}
                _source_.qualifiers = _source_QUALIFIERS
                _record_.features.append(_source_)

                """Line by line reading of the subsequence ids file, which contain the ids shared between the annotation
                and the gff file."""
                with open(SUB_SEQUENCES_FILE, "r") as sub_sequences:
                    sub_sequences_lines = (line.rstrip("\n") for line in sub_sequences.readlines())
                    try:
                        while True:
                            sub_sequence_id = next(sub_sequences_lines)

                            """Initialization of the minimized files for the features creation for each subsequence & features
                            initialization"""
                            SUB_SEQUENCE_FOLDER = SEQUENCE_FOLDER + sub_sequence_id + "/"
                            SUB_SEQUENCE_ANNOTATIONS_FILE = SUB_SEQUENCE_FOLDER + "anno.out"
                            SUB_SEQUENCE_GFF_FILE = SUB_SEQUENCE_FOLDER + "data.gff"
                            cmd("mkdir " + SUB_SEQUENCE_FOLDER)
                            cmd("grep \"" + sub_sequence_id + "\" " + SEQUENCE_ANNOTATIONS_FILE + " > " + SUB_SEQUENCE_ANNOTATIONS_FILE)
                            cmd("grep \"" + sub_sequence_id + "\" " + GFF_FILE + " > " + SUB_SEQUENCE_GFF_FILE)
                            assert os.path.isdir(SUB_SEQUENCE_FOLDER), "Could not create sub sequence folder!"
                            assert os.path.isfile(SUB_SEQUENCE_ANNOTATIONS_FILE), "Could not create sub sequence annotations file!"
                            assert os.path.isfile(SUB_SEQUENCE_GFF_FILE), "Could not create sub sequence gff file!"

                            _gene_LOCATION = None
                            _gene_QUALIFIERS = {"gene":list(), "note":list()}
                            _mRNA_LOCATION = list()
                            _mRNA_QUALIFIERS = {"gene":list(), "standard_name":list()}
                            _CDS_LOCATION = list()
                            _CDS_QUALIFIERS = {"gene":list(), "protein_id":list(), "product":list(), "note":list(), "db_xref":list(), "translation":list()}
                            _3_UTR_LOCATION = None
                            _5_UTR_LOCATION = None

                            """Features location creation with the gff file & Features qualifiers creation with the annotation file"""
                            with open(SUB_SEQUENCE_GFF_FILE, "r") as sub_sequence_gff:
                                sub_sequence_gff_lines = (line.rstrip("\n") for line in sub_sequence_gff.readlines())
                                try:
                                    while True:
                                        gff_entry = GFF_ENTRY(next(sub_sequence_gff_lines))

                                        if gff_entry.type == "mRNA": _gene_LOCATION = FeatureLocation(gff_entry.start, gff_entry.end, strand=gff_entry.strand)
                                        elif gff_entry.type == "exon": _mRNA_LOCATION.append(FeatureLocation(gff_entry.start, gff_entry.end, strand=gff_entry.strand))
                                        elif gff_entry.type == "CDS":  _CDS_LOCATION.append(FeatureLocation(gff_entry.start, gff_entry.end, strand=gff_entry.strand))
                                        elif gff_entry.type == "three_prime_UTR": _3_UTR_LOCATION = FeatureLocation(gff_entry.start, gff_entry.end, strand=gff_entry.strand)
                                        elif gff_entry.type == "five_prime_UTR":  _5_UTR_LOCATION = FeatureLocation(gff_entry.start, gff_entry.end, strand=gff_entry.strand)
                                except StopIteration:
                                    pass
                                finally:
                                    sub_sequence_gff.close()

                            with open(SUB_SEQUENCE_ANNOTATIONS_FILE, "r") as sub_sequence_annotations:
                                sub_sequence_annotations_lines = (line.rstrip("\n") for line in sub_sequence_annotations.readlines())
                                try:
                                    while True:
                                        annotation_entry = ANNOTATION_ENTRY(next(sub_sequence_annotations_lines))

                                        if annotation_entry.type == "qseq":  _CDS_QUALIFIERS["translation"].append(annotation_entry.desc)
                                        elif annotation_entry.type == "DE":
                                            _gene_QUALIFIERS["note"].append(annotation_entry.desc)
                                            _CDS_QUALIFIERS["product"].append(annotation_entry.desc)
                                        elif annotation_entry.type in ["BP_ARGOT", "CC_ARGOT", "MF_ARGOT"]: _CDS_QUALIFIERS["db_xref"].append("GO:" + annotation_entry.id)
                                        elif annotation_entry.type == "EC_ARGOT":  _CDS_QUALIFIERS["db_xref"] = _CDS_QUALIFIERS["db_xref"] + [x for x in annotation_entry.keggId()]
                                except StopIteration:
                                    pass
                                finally:
                                    sub_sequence_annotations.close()

                            """Features qualifiers appending and features appending to the record"""
                            _gene_QUALIFIERS["gene"] = list(set(_gene_QUALIFIERS["gene"]))
                            _gene_QUALIFIERS["note"] = list(set(_gene_QUALIFIERS["note"]))
                            _mRNA_QUALIFIERS["gene"] = list(set(_mRNA_QUALIFIERS["gene"]))
                            _mRNA_QUALIFIERS["standard_name"] = list(set(_mRNA_QUALIFIERS["standard_name"]))
                            _CDS_QUALIFIERS["gene"] = list(set(_CDS_QUALIFIERS["gene"]))
                            _CDS_QUALIFIERS["protein_id"] = list(set(_CDS_QUALIFIERS["protein_id"]))
                            _CDS_QUALIFIERS["product"] = list(set(_CDS_QUALIFIERS["product"]))
                            _CDS_QUALIFIERS["note"] = list(set(_CDS_QUALIFIERS["note"]))
                            _CDS_QUALIFIERS["db_xref"] = list(set(_CDS_QUALIFIERS["db_xref"]))
                            _CDS_QUALIFIERS["translation"] = list(set(_CDS_QUALIFIERS["translation"]))
                            _gene_ = SeqFeature(_gene_LOCATION, type="gene")
                            _gene_.qualifiers = _gene_QUALIFIERS
                            _mRNA_ = SeqFeature(mergeLocations(_mRNA_LOCATION), type="mRNA")
                            _mRNA_.qualifiers = _mRNA_QUALIFIERS
                            _CDS_ = SeqFeature(mergeLocations(_CDS_LOCATION), type="CDS")
                            _CDS_QUALIFIERS["db_xref"].sort()
                            _CDS_.qualifiers = _CDS_QUALIFIERS
                            _3_UTR_ = None if _3_UTR_LOCATION == None else SeqFeature(_3_UTR_LOCATION, type="3'UTR")
                            _5_UTR_ = None if _5_UTR_LOCATION == None else SeqFeature(_5_UTR_LOCATION, type="5'UTR")
                            _record_.features.append(_gene_)
                            _record_.features.append(_mRNA_)
                            _record_.features.append(_CDS_)
                            if _3_UTR_ != None: _record_.features.append(_3_UTR_)
                            if _5_UTR_ != None: _record_.features.append(_5_UTR_) 
                    except StopIteration:
                        pass
                    finally:
                        sub_sequences.close()

                """Exons features creation for the whole record"""
                with open(SEQUENCE_EXONS_FILE, "r") as exons:
                    exons_lines = (line.rstrip("\n") for line in exons.readlines())
                    try:
                        while True:
                            exon_entry = GFF_ENTRY(next(exons_lines))
                            #if exon_entry.length() < 15: continue
                            _exon_ = SeqFeature(FeatureLocation(exon_entry.start, exon_entry.end, exon_entry.strand), type="exon")
                            _exon_.qualifiers = {"note":list()}
                            _exon_.qualifiers["note"].append(exon_entry.exonId())
                            _record_.features.append(_exon_)
                    except StopIteration:
                        pass
                    finally:
                        exons.close()              
            #
                #SEQUENCE_MISC_FILE = SEQUENCE_FOLDER + "misc.gff"
                #cmd("grep \"" + sequence_id + "\" " + GFF_FILE + " | grep \":hit:\" > " + SEQUENCE_MISC_FILE)
                #assert os.path.isfile(SEQUENCE_MISC_FILE)
                #with open(SEQUENCE_MISC_FILE, "r") as misc_features:
                #    misc_features_lines = (line.rstrip("\n") for line in misc_features.readlines())
                #    try:
                #        while True:
                #            misc_feature_entry = GFF_ENTRY(next(misc_features_lines))
                #            _misc_feature_ = SeqFeature(FeatureLocation(misc_feature_entry.start, misc_feature_entry.end, misc_feature_entry.strand), type="misc_feature")
                #            _misc_feature_.qualifiers = {"note":list()}
                #            _record_.features.append(_misc_feature_)
                #    except StopIteration:
                #        pass
                #    finally:
                #        pass
                #    misc_features.close()

            except StopIteration:
                pass
            finally:
                cmd("rm -r " + SEQUENCE_FOLDER)

                """Outputing record in embl format """
                with open(OUTPUT_FILE, "w") as output_file:
                    print(_record_.format("embl"), file=output_file)

                cmd("java -jar embl-api-validator-1.1.265.jar " + OUTPUT_FILE)

        sequences.close()

if __name__ == "__main__":
    gffFile = fastaFile = annoFile = projFile = None

    args = sys.argv[1:]
    for i in [0,2,4,6]:
        if args[i] in ["-gff", "-g"]:
            gffFile = args[i+1]
        elif args[i] in ["-fasta", "-f"]:
            fastaFile = args[i+1]
        elif args[i] in ["-anno", "-a"]:
            annoFile = args[i+1]
        elif args[i] in ["-proj", "-p"]:
            projFile = args[i+1]
        else:
            assert False, "Invalid file input"

    main(fastaFile, annoFile, gffFile, projFile)