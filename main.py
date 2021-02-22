#!/usr/bin/python3
# -*-coding:utf-8

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

import os.path
import re
import subprocess
import sys

"""
The FASTA entries are created by the fastaIterator from fasta files.

Fasta entries are written as two lines handles with the sequence's id and the 
sequence what can be represented as follow
# >SEQ00001
# ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG...
"""
class FASTA_ENTRY():
    def __init__(self, id, seq):
        self.id = id
        self.seq = seq

"""
GFF entries are created by the gffIterator from gff files.

GFF entries are written as one line handles where fields are delimited by
tabulations as follow
# ptg001366l	.	contig	1	28502	.	.	.	ID=ptg01366l;Name=ptg01366l
      <0>      <1>   <2>   <3>   <4>   <5> <6> <7>              <8>

We select the required fields for the EMBL conversion, which are: 0, 2, 3, 4, 
5, 6, 8. The start position (field 3) is altered so that biopython don't ignore
the first position and the strand is converted in the +1/-1 format.

The exonId method allow to create an unique id for the exons which is required
for the EMBL conversion.
"""
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

"""
Annotation entries are created by the annotation itterator.

Annotation entries are written as one line handles where fields are delimited by 
tabulations as follow
# gene	EC_ARGOT	14.44	0.847	EC:4.2.1.134 	GO:0102345
  <0>   <1>         <2>     <3>     <4>             <5>   

We select the required fields for the EMBL conversion, which are 0, 1,4 and 5.

The keggId method allow the conversion of the various values values in EC_ARGOT
entry into one more relevant entry.
"""
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

"""
The four following generators are responsible to read and convert file entries
into their object type.
"""
def lineIterator(lineFile):
    assert os.path.isfile(lineFile)
    with open(lineFile, "r") as file:
        lines = (line.rsplit("\n")[0] for line in file.readlines())
        try:
            while True:
                (yield next(lines))
        except StopIteration:
            pass
        finally:
            file.close()

def fastaIterator(fastaFile):
    assert os.path.isfile(fastaFile)
    with open(fastaFile, "r") as file:
        lines = (line.rsplit("\n")[0] for line in file.readlines())
        try:
            while True:
                (yield FASTA_ENTRY(next(lines)[1:], next(lines)))
        except StopIteration:
            pass
        finally:
            file.close()

def gffIterator(gffFile):
    assert os.path.isfile(gffFile)
    with open(gffFile, "r") as file:
        lines = (line.rsplit("\n")[0] for line in file.readlines())
        try:
            while True:
                (yield GFF_ENTRY(next(lines)))
        except StopIteration:
            pass
        finally:
            file.close()

def annotationIterator(annotationFile):
    assert os.path.isfile(annotationFile)
    with open(annotationFile, "r") as file:
        lines = (line.rsplit("\n")[0] for line in file.readlines())
        try:
            while True:
                (yield ANNOTATION_ENTRY(next(lines)))
        except StopIteration:
            pass
        finally:
            file.close()

# Defining a "partial" method for subprocess.call
def cmd(_cmd_):
    subprocess.call(_cmd_, shell=True)

# Command responsible for the merging of Feature locations
def mergeLocations(_locationArray_):
    assert (len(_locationArray_) != 0), "Empty location array" 
    return _locationArray_[0] if len(_locationArray_) == 1 else CompoundLocation(_locationArray_)


def main(projectFile, fastaFile, gffFile, annotationFile):

    # Path to temp and output folder
    TEMP_FOLDER = "temp/"
    OUT_FOLDER = "out/"

    # Removing old temp & out folder if they already exists
    if os.path.isdir(TEMP_FOLDER): cmd("rm -r " + TEMP_FOLDER)
    if os.path.isdir(OUT_FOLDER): cmd("rm -r " + OUT_FOLDER)

    # Creating temp & out folder
    cmd("mkdir " + TEMP_FOLDER)
    cmd("mkdir " + OUT_FOLDER)

    # Path to submited files
    PROJECT_FILE = projectFile
    FASTA_FILE = fastaFile
    GFF_FILE = gffFile
    ANNOTATION_FILE = annotationFile

    # Verifying that folder are created and files exists
    assert os.path.isdir(TEMP_FOLDER), "Could not create temp folder!"
    assert os.path.isdir(OUT_FOLDER), "Could not create out folder!"
    assert os.path.isfile(FASTA_FILE), "FASTA file not found!"
    assert os.path.isfile(ANNOTATION_FILE), "Annotations file not found!"
    assert os.path.isfile(GFF_FILE), "GFF file not found!"
    assert os.path.isfile(projectFile), "Project file not found!"

    # Project variables initialization
    __DESCRIPTION__ = None
    __DIVISION__ = None
    __MOLECULE_TYPE__ = None
    __ORGANISM__ = None
    __PROJECT__ = None
    __TAXONOMY__ = None
    __TOPOLOGY__ = None

    for project_entry in (lineIterator(PROJECT_FILE)):
        parsed = re.split(r":", project_entry)
        if parsed[0] == "DESCRIPTION": __DESCRIPTION__ = parsed[1]
        elif parsed[0] == "DIVISION": __DIVISION__ = parsed[1]
        elif parsed[0] == "MOLECULE_TYPE": __MOLECULE_TYPE__ = parsed[1]
        elif parsed[0] == "ORGANISM": __ORGANISM__ = parsed[1]
        elif parsed[0] == "PROJECT": __PROJECT__ = parsed[1]
        elif parsed[0] == "TAXONOMY": __TAXONOMY__ = re.split(r"-",parsed[1])
        elif parsed[0] == "TOPOLOGY": __TOPOLOGY__ = parsed[1]

    # Verifying that the project variables are all defined
    assert __DESCRIPTION__, "Description undefined"
    assert __DIVISION__, "Division undefined"
    assert __MOLECULE_TYPE__, "Molecule type undefined"
    assert __ORGANISM__, "Organism undefined"
    assert __PROJECT__, "Project undefined"
    assert __TAXONOMY__ , "Taxonomy undefined"
    assert __TOPOLOGY__, "Topology undefined"

    # Contig file creation (the data contained in this file is complementary
    # to the data in the fasta file)
    CONTIG_FILE = TEMP_FOLDER + "contig.gff"
    cmd("grep \"contig\" " + GFF_FILE + " > " + CONTIG_FILE)
    assert os.path.isfile(CONTIG_FILE)

    # Iteration through the sequences from the fasta file
    for fasta_entry in fastaIterator(FASTA_FILE):

        # Subfiles and subfolder creation
        SEQUENCE_FOLDER = TEMP_FOLDER + fasta_entry.id + "/"
        cmd("mkdir " + SEQUENCE_FOLDER)
        assert os.path.isdir(SEQUENCE_FOLDER)

        SEQUENCE_GFF_FILE = SEQUENCE_FOLDER + "data.gff"
        SEQUENCE_EXON_FILE = SEQUENCE_FOLDER + "exons.gff"
        SEQUENCE_ANNOTATION_FILE = SEQUENCE_FOLDER + "anno.out"
        SEQUENCE_SUBSEQUENCE_FILE = SEQUENCE_FOLDER + "subseq.out"

        cmd("grep " + fasta_entry.id + " " + GFF_FILE + " > " + SEQUENCE_GFF_FILE)
        cmd("grep " + fasta_entry.id + " " + GFF_FILE + " | grep \":exon:\" > " + SEQUENCE_EXON_FILE)
        cmd("grep " + fasta_entry.id + " " + ANNOTATION_FILE + " > " + SEQUENCE_ANNOTATION_FILE)
        cmd("grep " + fasta_entry.id + " " + ANNOTATION_FILE + " | grep \"original_DE\" | cut -f1 > " + SEQUENCE_SUBSEQUENCE_FILE)

        # Assertion of the creation of the files
        assert os.path.isfile(SEQUENCE_GFF_FILE)
        assert os.path.isfile(SEQUENCE_EXON_FILE)
        assert os.path.isfile(SEQUENCE_ANNOTATION_FILE)
        assert os.path.isfile(SEQUENCE_SUBSEQUENCE_FILE)

        # Creation of the record output file
        OUTPUT_FILE = OUT_FOLDER + fasta_entry.id + ".dat"
        cmd("touch " + OUTPUT_FILE)
        assert os.path.isfile(OUTPUT_FILE)

        # Record creation
        _record_ = SeqRecord(
            Seq(fasta_entry.seq), 
            id=fasta_entry.id, 
            dbxrefs=["Project:" + __PROJECT__],
            annotations={
                "data_file_division":__DIVISION__, 
                "molecule_type":__MOLECULE_TYPE__, 
                "organism":__ORGANISM__, 
                "taxonomy":__TAXONOMY__,
                "topology":__TOPOLOGY__},
            description=__DESCRIPTION__)

        # Source feature creation
        source_entry = GFF_ENTRY(subprocess.check_output("grep \"" + fasta_entry.id + "\" " + CONTIG_FILE, shell=True).decode("utf-8"))
        _source_ = SeqFeature(
            FeatureLocation(source_entry.start, source_entry.end, strand=source_entry.strand),
            type="source",
            qualifiers={
                "organism":__ORGANISM__,
                "mol_type":__MOLECULE_TYPE__,
                "db_xref":list()
            })
        _record_.features.append(_source_)

        for sub_sequence_id in lineIterator(SEQUENCE_SUBSEQUENCE_FILE):
            
            # Subfiles and subfolder creation
            SUBSEQUENCE_FOLDER = SEQUENCE_FOLDER + sub_sequence_id + "/"
            cmd("mkdir " + SUBSEQUENCE_FOLDER)
            assert os.path.isdir(SUBSEQUENCE_FOLDER)

            SUBSEQUENCE_GFF_FILE = SUBSEQUENCE_FOLDER + "data.gff"
            SUBSEQUENCE_ANNOTATION_FILE = SUBSEQUENCE_FOLDER + "anno.out"

            cmd("grep " + sub_sequence_id + " " + SEQUENCE_GFF_FILE + " > " + SUBSEQUENCE_GFF_FILE)
            cmd("grep " + sub_sequence_id + " " + SEQUENCE_ANNOTATION_FILE + " > " + SUBSEQUENCE_ANNOTATION_FILE)

            # Assertion of the creation of the previous files and folders
            assert os.path.isfile(SUBSEQUENCE_GFF_FILE)
            assert os.path.isfile(SUBSEQUENCE_ANNOTATION_FILE)

            # Subsequence features initialization
            _gene_ = {"location":None, "qualifiers":{"gene":sub_sequence_id, "note":list()}, "type":"gene"}
            _mRNA_ = {"location":list(), "qualifiers":{"gene":sub_sequence_id, }, "type":"mRNA"}
            _CDS_ = {"location":list(), "qualifiers":{"gene":sub_sequence_id, "product":list(), "note":list(), "db_xref":list(), "translation":list()}, "type":"CDS"}
            _3UTR_ = {"location":None, "qualifiers":{"gene":sub_sequence_id, }, "type":"3'UTR"}
            _5UTR_ = {"location":None, "qualifiers":{"gene":sub_sequence_id, }, "type":"5'UTR"}

            # Subsequence features location
            for gff_entry in gffIterator(SUBSEQUENCE_GFF_FILE):
                location = FeatureLocation(gff_entry.start, gff_entry.end, strand=gff_entry.strand)
                if gff_entry.type == "mRNA": _gene_["location"] = location
                elif gff_entry.type == "exon": _mRNA_["location"].append(location)
                elif gff_entry.type == "CDS": _CDS_["location"].append(location)
                elif gff_entry.type == "three_prime_UTR": _3UTR_["location"] = location
                elif gff_entry.type == "five_prime_UTR": _5UTR_["location"] = location
            
            # Merging features location to produce CompoundLocation
            _mRNA_["location"] = mergeLocations(_mRNA_["location"])
            _CDS_["location"] = mergeLocations(_CDS_["location"])

            # Subsequence features annotations
            for annotation_entry in annotationIterator(SUBSEQUENCE_ANNOTATION_FILE):
                if annotation_entry.type == "qseq":
                    _CDS_["qualifiers"]["translation"].append(annotation_entry.desc)
                elif annotation_entry.type == "DE":
                    _gene_["qualifiers"]["note"].append(annotation_entry.desc)
                    _CDS_["qualifiers"]["product"].append(annotation_entry.desc)
                elif annotation_entry.type in ["BP_ARGOT", "CC_ARGOT", "MF_ARGOT"]:
                    _CDS_["qualifiers"]["db_xref"].append("GO:" + annotation_entry.id)
                elif annotation_entry.type == "EC_ARGOT":
                     _CDS_["qualifiers"]["db_xref"].extend(annotation_entry.keggId())

            # Appending features to record
            for feature in [_gene_, _mRNA_, _CDS_, _3UTR_, _5UTR_]:
                if feature["location"]: 
                    _record_.features.append(SeqFeature(feature["location"], type=feature["type"], qualifiers=feature["qualifiers"]))

        # Exon features creation
        for gff_entry in gffIterator(SEQUENCE_EXON_FILE):
            _record_.features.append(SeqFeature(FeatureLocation(gff_entry.start, gff_entry.end, strand=gff_entry.strand), type="exon", qualifiers={"note":[gff_entry.exonId()]}))

        # Printing output
        with open(OUTPUT_FILE, "w") as file:
            print(_record_.format("embl"), file=file)
            file.close()

        # Output validation
        cmd("java -jar embl-api-validator-1.1.265.jar " + OUTPUT_FILE)

        # Deleting temp sequence folder
        cmd("rm -r " + SEQUENCE_FOLDER)

        print("end of seq")


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

    main(projFile, fastaFile, gffFile, annoFile)   