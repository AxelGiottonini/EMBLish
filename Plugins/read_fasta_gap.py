#read_fasta_gap.py

import re

from Plugins.__read__ import __Read__
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

class Plugin(__Read__):
    
    """
    """
    def feature_initialize(self, pre_feature, metadata):
        return SeqFeature(
            FeatureLocation(pre_feature[0], pre_feature[1], 1),
            type="gap",
            qualifiers={})

    """
    """
    def multi_feature_initialize(self, pre_multi_feature, metadata):
        for element in pre_multi_feature:
            yield self.feature_initialize(element, metadata)
    
    """
    """
    def callbacks(self, *args):
        pass

    """
    """
    def process(self, app, caller_mode, key_handle, calls:list=[], target=None):
        try:
            temp = str(list(filter(lambda seq: seq.id == target, app.handles[key_handle]))[0].seq)
            pre_multi_feature = [gap.span() for gap in re.finditer(r"n+", temp)]
            
            feature = self.multi_feature_initialize(
                    pre_multi_feature,
                    app.metadata)
        except KeyError:
            return None

        receiver = self.callbacks(
                app,
                calls,
                target)

        return self.merge(feature, receiver)
