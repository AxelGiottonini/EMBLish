#verify_gff_maker_gene.py

from Plugins.__verify__ import __Verify__, FailedVerification
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

class Plugin(__Verify__):
    
    def process(self, app, element_to_verify):
        conversion_map = {"A":"T", "T":"A", "C":"G", "G":"C"}
        convert = lambda array : "".join([conversion_map[element] for element in array[::-1]])

        start_codon = None
        stop_codon = None

        feature_location = element_to_verify.location
        if(isinstance(feature_location, FeatureLocation)):
            if(feature_location.strand == 1):
                start_codon = app.current_sequence[feature_location._start:feature_location._start+3]
                stop_codon = app.current_sequence[feature_location._end-3:feature_location._end]
            else:
                start_codon = convert(app.current_sequence[feature_location._end-3:feature_location._end])
                stop_codon = convert(app.current_sequence[feature_location._start:feature_location._start+3])
        elif(isinstance(feature_location, CompoundLocation)):
            if(feature_location.strand == 1):
                start_codon = app.current_sequence[feature_location.parts[0]._start:feature_location.parts[0]._start+3]
                stop_codon = app.current_sequence[feature_location.parts[-1]._end-3:feature_location.parts[-1]._end]
            else:
                start_codon = convert(app.current_sequence[feature_location.parts[-1]._end-3:feature_location.parts[-1]._end])
                stop_codon = convert(app.current_sequence[feature_location.parts[0]._start:feature_location.parts[0]._start+3])

        if start_codon not in ["ATG"]:
            raise FailedVerification(f"invalid start codon: {start_codon}")

        if stop_codon not in ["TGA", "TAG", "TAA"]:
            raise FailedVerification(f"invalid stop codon: {stop_codon}")

        return None 