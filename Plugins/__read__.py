# __read__.py

from Plugins.__plugin__ import __Plugin__, RequiredMetadataError, UndefinedMethodError
from Plugins.__caller__ import Caller, CallerFailedVerification

class __Read__(__Plugin__):
    
    """
    """
    def feature_initialize(self, pre_feature, metadata):
        raise UndefinedMethodError("feature_initialize has not been defined.")

    """
    """
    def callbacks(self, app, calls, target):
        raise UndefinedMethodError("callbacks has not been defined.")

    """
    """
    def callbacks_extend(self, app, calls, target):
        caller = Caller(app)
        sender = []
        for app, key_plugin, *args in calls:
            temp = None
            try:
                temp = caller.run(app.plugins[key_plugin].process, app, *args, target)
            except CallerFailedVerification:
                sender = []
            if temp:
                sender.extend(temp)
        return sender

    """
    """
    def callbacks_append(self, app, calls, target):
        caller = Caller(app)
        sender = []
        for app, key_plugin, *args in calls:
            temp = None
            try:
                temp = caller.run(app.plugins[key_plugin].process, app, *args, target)
            except CallerFailedVerification:
                sender = []
            if temp:
                sender.append(temp)
        return sender

    """
    """
    def merge(self, feature, receiver):
        return feature

    """
    """
    def required_metadata_check(self, app, keys:list=[]):
        if keys:
            for key in keys:
                if not key in app.metadata:
                    raise RequiredMetadataError(f"Required metadata attribute, {key}, not found.")
        return True