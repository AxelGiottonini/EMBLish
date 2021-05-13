#read_tab_pannzer_CDS

import pandas as pd

class Plugin:

    def process(self, handle, metadata, calls:list=[], target=None):

        #initialisation
        try:
            anno_bp = handle.loc[(target[1], "BP_ARGOT"),:].reset_index()["id"]
        except KeyError:
            anno_bp = pd.Series([])

        try:
            anno_cc = handle.loc[(target[1], "CC_ARGOT"),:].reset_index()["id"]
        except KeyError:
            anno_cc = pd.Series([])

        try:
            anno_mf = handle.loc[(target[1], "MF_ARGOT"),:].reset_index()["id"]
        except KeyError:
            anno_mf = pd.Series([])

        try:
            anno_qsec = [handle.loc[(target[1], "qseq"),:].reset_index().iloc[0,1]]
        except KeyError:
            anno_qsec = list()

        try:
            anno_de =  [handle.loc[(target[1], "DE"),:].reset_index().iloc[0,1]]
        except KeyError:
            anno_de = list()

        _annotations_ = [{
            "db_xref":[f"GO:{str(go)}" for go in pd.concat([anno_bp, anno_cc, anno_mf])],
            "translation": anno_qsec,
            "product": anno_de
        }]

    
        #calls
        receiver = []
        for call,*args in calls:
            receiver.extend(call.process(*args, target=target))

        #output
        return _annotations_