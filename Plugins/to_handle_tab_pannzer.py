#tab_pannzer2handle.py

import pandas as pd
from Plugins.__plugin__ import __Plugin__

class Plugin(__Plugin__):

    def process(self, file_path):
        with open(file_path) as handle:
            anno = pd.read_csv(handle, sep="\t")

            anno.sort_values(by=["qpid", "type"], inplace=True)
            anno.drop(["score", "PPV"], axis=1, inplace=True)
            anno.set_index(["qpid", "type", "id", "desc"], inplace=True)
            
            return anno