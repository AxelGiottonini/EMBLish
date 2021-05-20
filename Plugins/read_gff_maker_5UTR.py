#read_gff_maker_5UTR.py

import pandas as pd
from Plugins.__read_gff_maker__ import __ReadGFFMaker__
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

class Plugin(__ReadGFFMaker__):

    """
    """
    def feature_initialize(self, pre_feature, metadata):
        return SeqFeature(
            FeatureLocation(int(pre_feature["start"]), int(pre_feature["stop"]), (1,-1)[pre_feature["strand"] == "-"]),
            type="5'UTR",
            qualifiers={
                "gene":None,
                "note":list()})
    
    """
    """
    def process(self, app, caller_mode, key_handle, calls:list=[], target=None):
        try:
            feature = self.feature_initialize(
                app.handles[key_handle].loc[(target[0], f"{target[1]}-mRNA-1", "five_prime_UTR"),:].reset_index().iloc[0,:],
                app.metadata)
        except KeyError:
            return None
            
        feature.qualifiers["gene"] = target[1]
        receiver = self.callbacks(
            app,
            calls,
            (target[0], f"{target[1]}-mRNA-1", "5'UTR"))
        
        return self.merge(feature, receiver)