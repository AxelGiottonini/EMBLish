#Plugins/read_gff_maker_exon.py

import pandas as pd

from Plugins.__read_gff_maker__ import __ReadGFFMaker__
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation


class Plugin(__ReadGFFMaker__):

    """
    """
    def feature_initialize(self, pre_feature, metadata):
        return SeqFeature(
            FeatureLocation(int(pre_feature["start"]), 
                            int(pre_feature["end"]), 
                            (1,-1)[pre_feature["strand"] == "-"]),
            type="exon",
            qualifiers={})
    
    """
    """
    def process(self, app, caller_mode, key_handle, calls:list=[], target=None):
        try:
            multi_pre_feature = app.handles[key_handle].loc[(target, slice(None), "exon"),:].iterrows()
        except KeyError:
            multi_pre_feature = []

        for _, pre_feature in multi_pre_feature:
            feature = self.feature_initialize(pre_feature, app.metadata)

            try:
                receiver = self.callbacks(
                        app,
                        calls,
                        target)
            except KeyError:
                receiver = []

            yield self.merge(feature, receiver)
