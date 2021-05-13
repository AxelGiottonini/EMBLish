#read_gff_maker_main.py

import pandas as pd
import re

from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

class Plugin:

    def process(self, handle, metadata, calls:list=[], target=None):
        
        for gene in handle.loc[(target, slice(None), "gene"),:].reset_index()["sub_seq_id"]:
            
            #initialize features
            _features_subset_ = []
               
            #calls
            receiver = []
            for call,*args in calls:
                receiver.extend(call.process(*args, target=(target, gene)))

            _features_subset_ = receiver
            yield _features_subset_