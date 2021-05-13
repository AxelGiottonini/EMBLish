#read_tab_pannzer_gene.py

#read_tab_pannzer_CDS

import pandas as pd

class Plugin:

    def process(self, handle, metadata, calls:list=[], target=None):

        #initialisation
        try:
            anno_de =  handle.loc[(target[1], "DE"),:].reset_index().iloc[0,1]
        except KeyError:
            anno_de = []

        _annotations_ = [{
            "note": anno_de
        }]

    
        #calls
        receiver = []
        for call,*args in calls:
            receiver.extend(call.process(*args, target=target))

        #output
        return _annotations_