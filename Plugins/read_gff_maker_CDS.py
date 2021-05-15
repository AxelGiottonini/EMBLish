#read_gff_maker_CDS.py

import pandas as pd
import itertools
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

class Plugin:

    """
    """
    def feature_initialize(self, pre_feature, metadata):
        refactor_pre_feature = lambda element: FeatureLocation(
            int(element[0]),
            int(element[1]),
            (1,-1)[element[2] == "-"])
        merge_pre_feature = lambda  array: array[0] if len(array) == 1 else CompoundLocation(array)

        return SeqFeature(
            merge_pre_feature(
                pre_feature.apply(refactor_pre_feature, axis=1)),
            type="CDS",
            qualifiers={
                "gene":None,
                "product":list(),
                "note":list(),
                "db_xref":list(),
                "translation":list(),
                "transl_table":metadata["transl_table"]})
    
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
        try:
            feature = self.feature_initialize(
                app.handles[key_handle].loc[(target[0], f"{target[1]}-mRNA-1", "CDS"),:].reset_index(),
                app.metadata)
        except KeyError:
            return None
            
        feature.qualifiers["gene"] = target[1]
        receiver = self.callbacks(
            app,
            calls,
            (target[0], f"{target[1]}-mRNA-1", "CDS"))
        
        return self.merge(feature, receiver)