#__plugin__.py

class __Plugin__:
    def process(*args):
        raise UndefinedMethodError("process has not been defined")

    def required_metadata_check(*args):
        return True

"""
"""
class RequiredMetadataError(Exception):
    def __init__(self, *args):
        if args:
            self.message = args[0]
        else:
            self.message = None

    def __str__(self):
        if self.message:
            return f"RequiredMetadataError: {self.message}"
        else:
            return "RequiredMetadataError has been raised"

"""
"""
class UndefinedMethodError(Exception):
    def __init__(self, *args):
        if args:
            self.message = args[0]
        else:
            self.message = None

    def __str__(self):
        if self.message:
            return f"UndefinedMethodError: {self.message}"
        else:
            return "UndefinedMethodError has been raised"