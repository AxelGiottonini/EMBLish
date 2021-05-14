#main.py

from core import app

from Bio import SeqIO
import pandas as pd
import re
import importlib

if __name__ == "__main__":

    _GLOBALS_ = {
        "handles":dict(),
        "plugins":dict(),
        "metadata":dict()
    }

    _PROCESSING_ = {
        "metadata":[],
        "plugins":[],
        "handles":[],
        "workflow":[]
    }
    current_field = None
    with open("files/config.info") as handle:
        for line in handle:
            current_line = line.rstrip("\n")

            if line[0] == "#": continue
            if current_line == "": continue

            if current_field:
                assert current_field in ["metadata", "plugins", "handles", "workflow"]
                if re.match(r"^(<\/)(\w+)(>)$", current_line):
                    current_field = None
                else:
                    _PROCESSING_[current_field].append(current_line) 
            else:
                assert re.match(r"^(<)(\w+)(>)$", current_line)
                assert current_line[1:-1] in ["metadata", "plugins", "handles", "workflow"]
                current_field = current_line[1:-1]
        handle.close()

    for element in _PROCESSING_["metadata"]:
            _GLOBALS_["metadata"][element.split(":")[0]] = element.split(":")[1]
    for element in _PROCESSING_["plugins"]:
            _GLOBALS_["plugins"][element.split(":")[0]] = importlib.import_module(element.split(":")[1].split(",")[0],element.split(":")[1].split(",")[1]).Plugin()
    for element in _PROCESSING_["handles"]:
            _GLOBALS_["handles"][element.split(":")[0]] = _GLOBALS_["plugins"][element.split(":")[1].split(",")[0]].process(element.split(":")[1].split(",")[1])

    max_level = 0
    for i in range(len(_PROCESSING_["workflow"])):
        element = _PROCESSING_["workflow"][i]
        regex = re.compile(r"^(-)+")
        level = len(regex.search(element).group())
        max_level = max(max_level, level)
        _PROCESSING_["workflow"][i] = (level, element[level:], [])

    _PROCESSING_["workflow"].insert(0, (0,None,[]))

    for i in range(max_level, -1, -1):
        for j in range(len(_PROCESSING_["workflow"])):
            element = _PROCESSING_["workflow"][j]
            if element[0] == i:
                for k in range(1, j+1):
                    if _PROCESSING_["workflow"][j-k] is not None and _PROCESSING_["workflow"][j-k][0] == i-1:
                        _PROCESSING_["workflow"][j-k][2].append((
                            _GLOBALS_["plugins"][_PROCESSING_["workflow"][j][1].split(",")[0]],
                            _GLOBALS_["handles"][_PROCESSING_["workflow"][j][1].split(",")[1]],
                            _GLOBALS_["metadata"],
                            _PROCESSING_["workflow"][j][2]
                        ))
                        _PROCESSING_["workflow"][j] = None
                        break
        _PROCESSING_["workflow"] = [x for x in _PROCESSING_["workflow"] if x is not None]

    app = app(_PROCESSING_["workflow"][0][2])
    app.run()