#Plugins/read_gff_maker_mRNA.py

import pandas as pd
from Plugins.__read_gff_maker__ import __ReadGFFMaker__
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

class Plugin(__ReadGFFMaker__):

    """
    """
    def feature_initialize(self, pre_feature, metadata):
        refactor_pre_feature = lambda element: FeatureLocation(
            int(element["start"]),
            int(element["end"]),
            (1,-1)[element["strand"] == "-"])
        merge_pre_feature = lambda  array: array[0] if len(array) == 1 else CompoundLocation(array)

        return SeqFeature(
            merge_pre_feature(pre_feature.apply(refactor_pre_feature, axis=1)),
            type="mRNA",
            qualifiers={
                "gene":None})

    """
    """
    def process(self, app, caller_mode, key_handle, calls:list=[], target=None):
        try:
            pre_feature = app.handles[key_handle].loc[(target[0], f"{target[1]}-mRNA-1:cds", "CDS"),:]
            feature = self.feature_initialize(
                    pre_feature,
                    app.metadata)
            feature.qualifiers["gene"] = target[1]
        except KeyError:
            return None

        try:
            receiver = self.callbacks(
                    app,
                    calls,
                    (target[0], f"{target[1]}-mRNA-1", "mRNA"))
        except KeyError:
            receiver = []

        return self.merge(feature, receiver)
