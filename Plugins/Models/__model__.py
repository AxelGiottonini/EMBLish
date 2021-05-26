#Plugins/Models/__model__.py

class __Model__():

    def __init__():
        raise UndefinedModelMethodError()


class UndefinedModelMethodError(Exception):
    def __init__(self, *args):
        if args:
            self.message = args[0]
        else self.message = None

    def __str__(self):
        if self.message:
            return f"UndefinedModelMethodError, {self.message}"
        else:
            return "UndefinedModelMethodError has been raised"
