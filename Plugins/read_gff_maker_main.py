#read_gff_maker_main.py

import pandas as pd

from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

class Plugin:
    """
    """
    def feature_initialize(self, pre_feature, metadata):
        return SeqFeature(
            FeatureLocation(int(pre_feature["start"]), int(pre_feature["stop"]), (1,-1)[pre_feature["strand"] == "-"]),
            type="source",
            qualifiers={
                "oganism":metadata["organism"],
                "mol_type":metadata["molecule_type"],
                "db_xref":list()})
    
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
    def callbacks_with_iterator(self, app, calls, target, iterator):
        sender = []

        for element in iterator:
            temp = self.callbacks(app, calls, (target, element))
            if temp:
                sender.extend(temp)
        return sender
            
    """
    """
    def merge(self, feature, receiver):
        return [feature] + receiver

    """
    """
    def process(self, app, key_handle, calls:list=[], target=None):
        try:
            feature = self.feature_initialize(
                app.handles[key_handle].loc[(target, slice(None), "contig"),:].reset_index().iloc[0,:], 
                app.metadata)
        except KeyError:
            return None

        try:
            receiver = self.callbacks_with_iterator(
                app, 
                calls, 
                target, 
                app.handles[key_handle].loc[(target, slice(None), "gene"),:].reset_index()["sub_seq_id"])
        except KeyError:
            receiver = []
        
        return  self.merge(feature, receiver)