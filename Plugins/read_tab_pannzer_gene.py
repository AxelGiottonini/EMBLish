#read_tab_pannzer_gene.py

import pandas as pd
from Plugins.__read_tab_pannzer__ import __ReadTabPannzer__

class Plugin(__ReadTabPannzer__):

    def feature_initialize(self, pre_feature, metadata):
        return {
            "note": self.feature_initialize_note(pre_feature)
        }

    def feature_initialize_note(self, pre_feature):
        sender = list()
        try:
            sender = [pre_feature("DE").iloc[0,1]]
        except KeyError:
            pass
        return sender

    def callbacks(self, app, calls, target):
        sender = []

        for app, key_plugin, *args in calls:
            temp = app.plugins[key_plugin].process(app, *args, target)
            if temp:
                sender += temp
    
        return sender

    def process(self, app, caller_mode, key_handle, calls:list=[], target=None):

        feature = self.feature_initialize(            
            (lambda field: app.handles[key_handle].loc[(target[1], field)].reset_index()),
            app.metadata)
        receiver = self.callbacks(
            app,
            calls,
            target
        )

        return self.merge(feature, receiver)