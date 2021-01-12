# Conversion and merging of MAKER2 output with PANNZER output

## Background

Genome annotation consists in two consecutive steps. We initially predict gene position and the we add functional annotations. Two files are obtained from these steps, one in the GFF3 format and one in a tabular format.

As those files are not suitable to be submitted to a database, we have to convert them in a standardized format. We thus decided to merge the results in the EMBL format, what is not restrictive as many converter exists based on this format.

One may also notice that similar converters already exists but they only aim to convert the GFF3 output in the EMBL format. Two approaches to our script are thus possible:

1. Converting the tabular format into the GFF3 format and then using an already existing API to convert both GFF3 files into an EMBL data entry.
2. Reading GFF3 and tabular files in parallel and creating an EMBL entry from the reads.

We decided to create our script with the second approach.

## Tools

The creation of the EMBL entry is done with the BioPython package, especially the SeqRecord class to create the entry and the SeqFeature class to describe in details our entry. An example of the code is shown below:

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

## Creating sequence
sequence = Seq("#")

## Creating record
record = SeqRecord(sequence, 
                   id="#")
record.dbxrefs = ["Project:#"]
record.annotations = {
    "data_file_division":"#"
    "molecule_type":"#"
    "organism":"#"
    "taxonomy":"#"
    "topology":"#"
}

## Creating features
featureA = SeqFeature(
    FeatureLocation(start, end, strand), 
    type="#")
...
featureZ = SeqFeature(...)

## Appending features to record
record.features.append(featureA)
...
record.features.append(featureZ)
```

In order to parse the data, we mainly use the regular expression package from python (`re`) and the `grep` function from Unix, executed through the `subprocess.call()` method.

Finally, aiming to validate the EMBL entry, we call the `embl-api-validator` from the ENA.

## Work flow

### Identifying the sequences

The sequence identification is done through the reading of the submitted FASTA file containing the sequences. The format of an entry is the following:

```
>seq_id
ATCGATCGATCGAT...
```

We thus have to read two lines at a time in order to obtain the sequence and its identifier. Once these data are obtained, we can create the `SeqRecord` object and carry the work by adding features to the record.

### Identifying the subsequences

In order to obtain the identifier of the subsequences, we filter the annotation file twice. The first filtering allows us the select the entries about the current sequence when the second filtering allows us to select unique subsequence identifiers. The result of the filtering is a file formatted as follow:

```
subsequence_id_A
...
subsequence_id_Z
```

Once we have got this file, we can start creating various types of features: `gene`, `mRNA`, `CDS`, `3'UTR`, `5'UTR`. The location of each feature is retrieved from the GFF3 file when the annotations are researched on the annotation file.

### Identifying the exons

The creation of the exon features is analogous to the subsequences workflow with also a two step filtering, one for the sequence, the second for the exon keyword but this time directly on the GFF3 file, and a file output. The outputted file is read line by line and each exon entry is converted into a record feature.

### Summary

The work flow is nested in loops, which may be summarized as follow:

```pseudocode
foreach sequence in fasta file:
	new record
	foreach subsequence in subsequences file
		new feature
		feature set location from GFF3
		feature set annotations from annotations file
		record append feature
	foreach exon in exons file
		new feature
		feature set location from GFF3
		feature set qualifiers from GFF3
		record append feature
```

## Feature definitions

In our script we had to define seven features represented here as `JSON` objects:

```json
//source
{
    "location":FeatureLocation(),
    "qualifiers":{
    	"organism":"#",
    	"mol_type":"#",
        "db_xref":[]
	}
}
//gene
{
    "location":FeatureLocation(),
    "qualifiers":{
        "gene":"#",
        "note":"#",
    }
}
//mRNA
{
    "location":[FeatureLocation()],
    "qualifiers":{
        "gene":"#",
        "standard_name":"#"
    }
}
//CDS
{
    "location":[FeatureLocation()],
    "qualifiers":{
        "gene":"#",
        "protein_id":"#",
        "note":"#",
        "db_xref":[],
        "translation":"#"
    }
}
//3'UTR & 5'UTR
{
    "location":FeatureLocation()
}
//exon
{
    "location":FeatureLocation(),
    "qualifiers":{
    	"note":"#"
	}
}
```

