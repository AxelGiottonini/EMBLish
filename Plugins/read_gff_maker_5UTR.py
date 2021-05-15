#read_gff_maker_5UTR.py

import pandas as pd
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

class Plugin:

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
    def callbacks(self, app, calls, target):
        sender = []

        for app, key_plugin, *args in calls:
            temp = app.plugins[key_plugin].process(app, *args, target)
            if temp:
                sender.append(temp)

        return sender

    """
    """
    def merge(self, feature, receiver):
        return feature

    """
    """
    def process(self, app, key_handle, calls:list=[], target=None):
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