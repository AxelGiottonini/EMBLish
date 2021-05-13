#read_gff_maker_3UTR.py

import pandas as pd
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

class Plugin:

    def process(self, handle, metadata, calls:list=[], target=None):
        try:
            location = (handle.loc[(target[0], f"{target[1]}-mRNA-1", "three_prime_UTR"),:].reset_index())

            _sub_features_ = [
                SeqFeature(
                    FeatureLocation(int(location.iloc[0,0]), int(location.iloc[0,1]), (1,-1)[location.iloc[0,2] == "-"]),
                    type="3'UTR",
                    qualifiers={
                        "gene":target[1],
                        "note":list()})]
                
            #calls
            receiver = []
            for call,*args in calls:
                receiver.extend(call.process(*args, target=(target[0], f"{target[1]}-mRNA-1", "3'UTR")))
            
            return _sub_features_
        except KeyError:
            return []