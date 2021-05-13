#read_gff_maker_gene.py

import pandas as pd
import itertools
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

class Plugin:

    def process(self, handle, metadata, calls:list=[], target=None):
        location = (handle.loc[(target[0], target[1], "gene"),:].reset_index())

        _sub_features_ = [
            SeqFeature(
                FeatureLocation(int(location.iloc[0,0]), int(location.iloc[0,1]), (1,-1)[location.iloc[0,2] == "-"]),
                type="gene",
                qualifiers={
                    "gene":target[1],
                    "note":list()})]
            
        #calls
        receiver = []
        for call,*args in calls:
            receiver.extend(call.process(*args, target=(target[0], target[1], "gene")))
        
        annotations = list(itertools.chain(receiver))
        for annotation in annotations:
            if "note" in annotation.keys() and annotation["note"] != []:
                for sub_feature in _sub_features_:
                    sub_feature.qualifiers["note"].extend(annotation["note"])

        return _sub_features_