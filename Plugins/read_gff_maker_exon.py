#read_gff_maker_exon.py

import pandas as pd
import re

from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation


class Plugin:

    def process(self, handle, metadata, calls:list=[], target=None):
        exons = handle.loc[(target, slice(None), "exon"),:].reset_index()

        for index, exon in exons.iterrows():
            _features_subset_ = [
                SeqFeature(
                    FeatureLocation(int(exon["start"]), int(exon["stop"]), (1,-1)[exon["strand"] == "-"]),
                    type="exon",
                    qualifiers={}
                )]

            #calls
            receiver = []
            for call, *args in calls:
                receiver.extend(call.process(*args, target=(target)))

             
            yield _features_subset_
