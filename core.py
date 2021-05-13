#core.py

import importlib

class app:
    def __init__(self, plugins:list=[]):
        assert plugins != [], "No plugins specified"

        self.plugins = plugins

    def run(self):
        for plugin,*args in self.plugins:
            plugin.process(*args)