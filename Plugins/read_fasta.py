#plugin.py

import itertools

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord 

class Plugin:

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
        sender = []

        for app, key_plugin, *args in calls:
            temp = app.plugins[key_plugin].process(app, *args, target)
            if temp:
                sender += temp
    
        return sender

    """
    """
    def merge(self, feature, receiver):
        feature.features = receiver    

    """
    """
    def process(self, app, key_handle, calls:list=[], target=None):

        for element in app.handles[key_handle]:
            feature = self.feature_initialize(element, app.metadata)
            receiver = self.callbacks(app, calls, (feature.id))
            self.merge(feature, receiver)
            
            with open(f"out/{feature.id}.dat", "w") as o:
                print(feature.format("embl"), file=o)