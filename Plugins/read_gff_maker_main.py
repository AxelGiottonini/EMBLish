#read_gff_maker_main.py

import pandas as pd
import re

from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

class Plugin:

    def process(self, handle, metadata, calls:list=[], target=None):
        
        location = (handle.loc[(target, slice(None), "contig"),:].reset_index())
        _feature_ = [
            SeqFeature(
                FeatureLocation(int(location.iloc[0,3]), int(location.iloc[0,4]), (1,-1)[location.iloc[0,5] == "-"]),
                type="source",
                qualifiers={
                    "oganism":metadata["organism"],
                    "mol_type":metadata["molecule_type"],
                    "db_xref":list()})]

        yield _feature_

        for gene in handle.loc[(target, slice(None), "gene"),:].reset_index()["sub_seq_id"]:
            
            #initialize features
            _features_subset_ = []
               
            #calls
            receiver = []
            for call,*args in calls:
                receiver.extend(call.process(*args, target=(target, gene)))

            _features_subset_ = receiver
            yield _features_subset_