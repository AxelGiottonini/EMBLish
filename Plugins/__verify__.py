#__verify__.py

from Plugins.__plugin__ import __Plugin__, RequiredMetadataError, UndefinedMethodError

class __Verify__(__Plugin__):
    pass

class FailedVerification(Exception):
    def __init__(self, *args):
        if args:
            self.message = args[0]
        else:
            self.message = None

    def __str__(self):
        if self.message:
            return f"FailedVerification: {self.message}"
        else:
            return "FailedVerification has been raised"