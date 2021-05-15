#read_gff_maker_exon.py

import pandas as pd
import re

from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation


class Plugin:

    """
    """
    def feature_initialize(self, pre_feature, metadata):
        return SeqFeature(
            FeatureLocation(int(pre_feature["start"]), int(pre_feature["stop"]), (1,-1)[pre_feature["strand"] == "-"]),
            type="exon",
            qualifiers={})
    
    """
    """
    def multi_feature_initialize(self, pre_multi_feature, metadata):
        for _, element in pre_multi_feature:
            yield self.feature_initialize(element, metadata)

    """
    """
    def callbacks(self, app, calls, target):
        sender = []

        for app, key_plugin, *args in calls:
            temp = app.plugins[key_plugin].process(app, *args, target)
            if temp:
                sender.append()

        return sender
    
    """
    """
    def merge(self, feature, receiver):
        return feature

    """
    """
    def process(self, app, key_handle, calls:list=[], target=None):
        try:
            feature = self.multi_feature_initialize(
                app.handles[key_handle].loc[(target, slice(None), "exon"),:].reset_index().iterrows(),
                app.metadata)
        except KeyError:
            return None
            
        receiver = self.callbacks(
            app,
            calls,
            target)
        return self.merge(feature, receiver)
