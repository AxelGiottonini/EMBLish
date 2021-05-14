#tab_pannzer2handle.py

import pandas as pd

class Plugin:

    def process(self, file_path):
        with open(file_path) as handle:
            anno = pd.read_csv(handle, sep="\t")
            anno = anno.sort_values(by=["qpid", "type"]).drop(["score", "PPV"], axis=1).set_index(["qpid", "type", "id", "desc"])
            return anno