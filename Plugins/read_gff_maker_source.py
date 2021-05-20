#read_gff_maker_source.py

import pandas as pd

from Plugins.__read_gff_maker__ import __ReadGFFMaker__
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

class Plugin(__ReadGFFMaker__):

    """
    """
    def feature_initialize(self, pre_feature, metadata):
        return SeqFeature(
            FeatureLocation(int(pre_feature["start"]), int(pre_feature["stop"]), (1,-1)[pre_feature["strand"] == "-"]),
            type="source",
            qualifiers={
                "organism":metadata["organism"],
                "mol_type":metadata["molecule_type"],
                "db_xref":list()})

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
    def process(self, app, caller_mode, key_handle, calls:list=[], target=None):
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

    """
    """
    def required_metadata_check(self, app, keys:list=[]):
        return super().required_metadata_check(app, ["organism", "molecule_type"])