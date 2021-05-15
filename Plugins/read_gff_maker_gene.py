#read_gff_maker_gene.py

import pandas as pd
import itertools
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

class Plugin:

    """
    """
    def feature_initialize(self, pre_feature, metadata):
        return SeqFeature(
            FeatureLocation(int(pre_feature["start"]), int(pre_feature["stop"]), (1,-1)[pre_feature["strand"] == "-"]),
            type="gene",
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
        for element in receiver:
            for key in element.keys():
                feature.qualifiers[key].extend(element[key])
        return feature

    """
    """
    def process(self, app, key_handle, calls:list=[], target=None):
        feature = self.feature_initialize(
            app.handles[key_handle].loc[(target[0], target[1], "gene"),:].reset_index().iloc[0,:],
            app.metadata)
        feature.qualifiers["gene"]=target[1]
        receiver = self.callbacks(
            app,
            calls,
            (target[0], f"{target[1]}-mRNA-1", "gene"))

        return self.merge(feature, receiver)