#plugin.py

import itertools

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord 
from Plugins.__read__ import __Read__

class Plugin(__Read__):
    """
    """
    def feature_initialize(self, pre_feature, metadata):
        return SeqRecord(
            pre_feature.seq,
            pre_feature.id,
            dbxrefs=["Project:" + metadata["project"]],
            annotations={
                "division":metadata["division"],
                "molecule_type":metadata["molecule_type"],
                "organism":metadata["organism"],
                "taxonomy":metadata["taxonomy"],
                "topology":metadata["topology"]},
            description="")
    
    """
    """
    def callbacks(self, app, calls, target):
        return super().callbacks_extend(app, calls, target)

    """
    """
    def merge(self, feature, receiver):
        feature.features = receiver    

    """
    """
    def process(self, app, caller_mode, key_handle, calls:list=[], target=None):

        for element in app.handles[key_handle]:
            feature = self.feature_initialize(element, app.metadata)
            app.current_sequence = feature.seq
            receiver = self.callbacks(app, calls, (feature.id))
            self.merge(feature, receiver)
            
            #with open(f"out/{feature.id}.dat", "w") as o:
            with open(f"{feature.id}.dat", "w") as o:
                print(feature.format("embl"), file=o)

    """
    """
    def required_metadata_check(self, app, keys:list=[]):
        return super().required_metadata_check(app, ["project", "transl_table", "molecule_type", "organism", "taxonomy", "topology"])
