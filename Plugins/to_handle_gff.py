#Plugins/to_handle_gff.py

import pandas as pd
import re
import sys
from functools import partial
from Plugins.__plugin__ import __Plugin__

class Plugin(__Plugin__):
    """
        Attributes
            file_path : string, path to the file to convert
        Return
            pandas data frame of 18 columns:
                 0. seqid
                 1. source
                 2. type
                 3. start
                 4. end
                 5. score
                 6. strand
                 7. phase
                 8. attributes_id
                 9. attributes_name
                10. attributes_alias
                11. attributes_parent
                12. attributes_target
                13. attributes_gap
                14. attributes_derives_from
                15. attributes_note
                16. attributes_dbxref
                17. attributes_ontology_term
    """
    def process(self, file_path):

        with open(file_path) as handle:
            temp = pd.read_csv(handle, sep="\t", comment="#", header=None)
            temp.columns = [
                "seqid",
                "source",
                "type",
                "start",
                "end",
                "score",
                "strand",
                "phase",
                "attributes"]
            
            temp.drop_duplicates(
                    subset=["seqid", "type", "start", "end", "strand"],
                    inplace=True)

            temp["start"] = temp["start"].apply(lambda x: x - 1)
           
            def extract_attribute(prefix, attributes):
                search = re.search(
                        rf"(?P<head>.*{prefix}=)+(?P<core>[^;]*)(?P<tail>;.*)*",
                        attributes)
                
                if search: return search.groupdict()["core"]
                
                return None

            temp["attributes_id"]               = temp["attributes"].apply(lambda x: extract_attribute("ID",x))
            temp["attributes_name"]             = temp["attributes"].apply(lambda x: extract_attribute("Name",x))
            temp["attributes_alias"]            = temp["attributes"].apply(lambda x: extract_attribute("Alias",x))
            temp["attributes_parent"]           = temp["attributes"].apply(lambda x: extract_attribute("Parent",x)) 
            temp["attributes_target"]           = temp["attributes"].apply(lambda x: extract_attribute("Target",x))
            temp["attributes_gap"]              = temp["attributes"].apply(lambda x: extract_attribute("Gap",x))
            temp["attributes_derives_from"]     = temp["attributes"].apply(lambda x: extract_attribute("Derives_from",x))
            temp["attributes_note"]             = temp["attributes"].apply(lambda x: extract_attribute("Note",x))
            temp["attributes_dbxref"]           = temp["attributes"].apply(lambda x: extract_attribute("Dbxref",x))
            temp["attributes_ontology_term"]    = temp["attributes"].apply(lambda x: extract_attribute("Ontology_term",x))

            temp.drop("attributes", axis=1, inplace=True)
            return temp
