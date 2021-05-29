#gff_maker2handle.py

import pandas as pd
import re
from Plugins.to_handle_gff import Plugin

class Plugin(Plugin):

    """
        Attributes
            file_path : string, path to the file to convert
        Return
            pandas data frame of __ columns:
                0. Index:
                    seqid
                    attributes_id
                    type
                1. start
                2. end
                3. strand
    """
    def process(self, file_path):
        
        temp = super().process(file_path)
        temp.drop(["source",
                   "score",
                   "phase",
                   "attributes_name",
                   "attributes_alias",
                   "attributes_parent",
                   "attributes_target",
                   "attributes_gap",
                   "attributes_derives_from",
                   "attributes_note",
                   "attributes_dbxref",
                   "attributes_ontology_term"], 
                   axis=1, 
                   inplace=True)
        
        temp = temp[["seqid", "attributes_id", "type", "start", "end", "strand"]]
        temp.sort_values(by=["seqid", "attributes_id", "type"],
                         inplace=True)
        temp.set_index(["seqid", "attributes_id", "type"],
                       inplace=True)
        return temp
