#read_gff_maker_CDS.py

import pandas as pd
import itertools
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

def mergeLocations(_locationArray_):
    return _locationArray_[0] if len(_locationArray_) == 1 else CompoundLocation(_locationArray_)

class Plugin:

    def process(self, handle, metadata, calls:list=[], target=None):
        locations = (handle.loc[(target[0], f"{target[1]}-mRNA-1", "CDS"),:].reset_index())

        _sub_features_ = [
            SeqFeature(
                mergeLocations(locations.apply(lambda location: FeatureLocation(int(location[0]), int(location[1]), (1,-1)[location[2] == "-"]), axis=1)),
                type="CDS",
                qualifiers={
                    "gene":target[1],
                    "product":list(),
                    "note":list(),
                    "db_xref":list(),
                    "translation":list(),
                    "transl_table":metadata["transl_table"]})]

        #calls
        receiver = []
        for call,*args in calls:
            receiver.extend(call.process(*args, target=(target[0], f"{target[1]}-mRNA-1", "CDS")))
        
        annotations = list(itertools.chain(receiver))
        for annotation in annotations:
            if "product" in annotation.keys() and annotation["product"] != []:
                for sub_feature in _sub_features_:
                    sub_feature.qualifiers["product"].extend(annotation["product"])

            if "note"  in annotation.keys() and annotation["note"] != []:
                for sub_feature in _sub_features_:
                    sub_feature.qualifiers["note"].extend(annotation["note"])

            if "db_xref" in annotation.keys() and annotation["db_xref"] != []:
                for sub_feature in _sub_features_:
                    sub_feature.qualifiers["db_xref"].extend(annotation["db_xref"])

            if "translation" in annotation.keys() and annotation["translation"] != []:
                for sub_feature in _sub_features_:
                    sub_feature.qualifiers["translation"].extend(annotation["translation"])

        return _sub_features_