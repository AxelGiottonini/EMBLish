#read_gff_maker_misc_feature.py

import pandas as pd
import itertools
from Plugins.__read_gff_maker__ import __ReadGFFMaker__
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

class Plugin(__ReadGFFMaker__):
    """
    """
    def feature_initialize(self, pre_feature, metadata):
        return SeqFeature(
            FeatureLocation(int(pre_feature["start"]), int(pre_feature["stop"]), (1,-1)[pre_feature["strand"] == "-"]),
            type="misc_feature",
            qualifiers={
                #"gene":None,
                "note":list()})

    """
    """
    def process(self, app, caller_mode, key_handle, calls:list=[], target=None):
        feature = self.feature_initialize(
            app.handles[key_handle].loc[(target[0], target[1], "gene"),:].reset_index().iloc[0,:],
            app.metadata)
        feature.qualifiers["note"].append(target[1])
        #feature.qualifiers["gene"]=target[1]
        receiver = self.callbacks(
            app,
            calls,
            (target[0], f"{target[1]}-mRNA-1", "misc_feature"))
            
        return self.merge(feature, receiver)