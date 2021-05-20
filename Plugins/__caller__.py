#__caller__.py

"""
"""
from Plugins.__verify__ import FailedVerification


class Caller:
    """
    """
    def __init__(self, app):
        self.app = app
        self.status = True
        self.last_result = None

    """
    """
    def run(self, callback_function, *args):
        if args[1] not in ["default", "verify", "bypass"]:
            raise UnknownCallerModeError(f"{args[1]} mode is not defined.")

        if self.status and args[1] in ["default"]:
            self.last_result = callback_function(*args)
            return self.last_result

        elif self.status and args[1] in ["verify"]:
            try:
                callback_function(self.app, self.last_result)
            except FailedVerification:
                self.status = False
                raise CallerFailedVerification()
            return None

        elif not self.status and args[1] in ["bypass"]:
            self.last_result = callback_function(*args)
            return self.last_result

class UnknownCallerModeError(Exception):
    def __init__(self, *args):
        if args:
            self.message = args[0]
        else:
            self.message = None

    def __str__(self):
        if self.message:
            return f"UnknownCallerModeError: {self.message}"
        else:
            return "UnknownCallerModeError has been raised"

class CallerFailedVerification(Exception):
    def __init__(self, *args):
        if args:
            self.message = args[0]
        else:
            self.message = None

    def __str__(self):
        if self.message:
            return f"CallerFailedVerification: {self.message}"
        else:
            return "CallerFailedVerification has been raised"