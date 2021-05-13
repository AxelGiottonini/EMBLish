#plugin.py

import importlib
import itertools

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord 

class Plugin:

    def process(self, handle, metadata, calls:list=[], target=None):

            for record in handle:

                #initialize record
                _record_ = SeqRecord(
                    record.seq,
                    record.id,
                    dbxrefs=["Project:" + metadata["project"]],
                    annotations={"division":metadata["division"],"molecule_type":metadata["molecule_type"],"organism":metadata["organism"],"taxonomy":metadata["taxonomy"],"topology":metadata["topology"]},
                    description=""
                )

                #calls
                receiver = []
                for call,*args in calls:
                    receiver.extend(call.process(*args, target=_record_.id))
        
                #post output treatment
                _record_.features = list(itertools.chain(*receiver))

                #outputing
                with open(f"out/{_record_.id}.dat", "w") as o:
                    print(_record_.format("embl"), file=o)