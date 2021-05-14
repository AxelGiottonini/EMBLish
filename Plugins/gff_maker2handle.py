#gff_maker2handle.py

import pandas as pd
import re

class Plugin:

    def process(self, file_path):
        with open(file_path) as handle:
            gff = pd.read_csv(handle, sep="\t")
            gff = gff.reset_index()
            gff.columns = ["seq_id", "source", "ft_type", "start", "stop", "score", "strand", "phase", "attr"]
            gff = gff.sort_values(by=["seq_id"]).drop(["source", "score", "phase"], axis=1).dropna()
            gff["sub_seq_id"] = [re.split(r':',re.search("^ID=.*?;", x).group(0)[3:-1])[0] for x in gff["attr"]]
            return gff[["seq_id", "sub_seq_id", "ft_type", "start", "stop", "strand"]].sort_values(by=["seq_id", "sub_seq_id", "ft_type"]).set_index(["seq_id", "sub_seq_id", "ft_type", "start", "stop", "strand"]) 
        