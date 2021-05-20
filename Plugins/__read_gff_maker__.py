#__read_gff_maker__.py

from Plugins.__read__ import __Read__

class __ReadGFFMaker__(__Read__):

    """
    """
    def multi_feature_initialize(self, pre_multi_feature, metadata):
        raise UndefinedMethodError("multi_feature_initialize has not been defined.")

    """
    """
    def callbacks(self, app, calls, target):
        return super().callbacks_append(app, calls, target)

    """
    """
    def callbacks_with_iterator(self, app, calls, target, iterator):
        raise UndefinedMethodError("callbacks_with_iterator has not been defined.")

