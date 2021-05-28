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
    def multi_feature_initialize(self, pre_multi_feature, metadata):
        for _, element in pre_multi_feature:
            yield self.feature_initialize(element, metadata)


    """
    """
    def process(self, app, caller_mode, key_handle, calls:list=[], target=None):
        try:
            pre_multi_feature = app.handles[key_handle].loc[(target, slice(None), "exon"),:].iterrows()
            feature = self.multi_feature_initialize(pre_multi_feature, app.metadata)
        except KeyError:
            return None
            
        receiver = self.callbacks(
            app,
            calls,
            target)
        return self.merge(feature, receiver)
