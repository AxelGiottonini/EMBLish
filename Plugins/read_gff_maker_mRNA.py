#read_gff_maker_mRNA.py

import pandas as pd
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

def mergeLocations(_locationArray_):
    return _locationArray_[0] if len(_locationArray_) == 1 else CompoundLocation(_locationArray_)

class Plugin:

    def process(self, handle, metadata, calls:list=[], target=None):
        locations = (handle.loc[(target[0], f"{target[1]}-mRNA-1", "CDS"),:].reset_index())

        #initialise
        _sub_features_ = [
            SeqFeature(
                mergeLocations(locations.apply(lambda location: FeatureLocation(int(location[0]), int(location[1]), (1,-1)[location[2] == "-"]), axis=1)),
                type="mRNA",
                qualifiers={
                    "gene":target[1]})]

        #calls
        receiver = []
        for call,*args in calls:
            receiver.extend(call.process(*args, target=(target[0], f"{target[1]}-mRNA-1", "mRNA")))
        
        return _sub_features_