# EMBLish

<i>
EMBLish aims to provide a tool for converting gene prediction files such as Maker output file and gene annotation files such as Pannzer2 output file into the EMBL format
 </i>
 
## Prerequistes:
- Python >= 3.9.0
- Biopython >= 1.78

## Download with git:
```
git clone https://github.com/AxelGiottonini/EMBLish.git
```

## Usage:

### Configuration file:

Configuration files are divided into four parts: `metadata`, `plugins`, `handles` and `workflow`. Each part is mandatory and we will present it below. User can also add comment lines using the `#` symbol.

The example we use here aims to convert a Maker output file wiht a Pannzer output file into an EMBL entry. The current directory contains four files : `config.info`, `sequences.fasta`, `data.gff`, `anno.out`.

#### `config.info`
```
# @author : Axel Giottonini
# @file : config.info
# @date : 19.05.2021

<metadata>
project:PRJEB1234
division:inv
taxonomy:29031
organism:Phlebotomus papatasi
molecule_type:genomic DNA
topology:linear
description:empty
transl_table:1,int
</metadata>

<plugins>
read_fasta:.read_fasta,Plugins
read_gff_maker_3UTR:.read_gff_maker_3UTR,Plugins
read_gff_maker_5UTR:.read_gff_maker_5UTR,Plugins
read_gff_maker_CDS:.read_gff_maker_CDS,Plugins
read_gff_maker_exon:.read_gff_maker_exon,Plugins
read_gff_maker_gene:.read_gff_maker_gene,Plugins
read_gff_maker_misc_feature:.read_gff_maker_misc_feature,Plugins
read_gff_maker_mRNA:.read_gff_maker_mRNA,Plugins
read_gff_maker_source:.read_gff_maker_source,Plugins
read_tab_pannzer_CDS:.read_tab_pannzer_CDS,Plugins
read_tab_pannzer_gene:.read_tab_pannzer_gene,Plugins
to_handle_fasta:.to_handle_fasta,Plugins
to_handle_gff_maker:.to_handle_gff_maker,Plugins
to_handle_tab_pannzer:.to_handle_tab_pannzer,Plugins
verify_gff_maker_CDS:.verify_gff_maker_CDS,Plugins
</plugins>

<handles>
fasta:to_handle_fasta,sequences.fasta
gff_maker:to_handle_gff_maker,data.gff
tab_pannzer:to_handle_tab_pannzer,anno.out
</handles>

<workflow>
-read_fasta,fasta
--read_gff_maker_source,gff_maker
---read_gff_maker_gene,gff_maker
----read_tab_pannzer_gene,tab_pannzer
---read_gff_maker_mRNA,gff_maker
---read_gff_maker_CDS,gff_maker
----read_tab_pannzer_CDS,tab_pannzer
---verify_gff_maker_CDS,no_file,verify
---read_gff_maker_misc_feature,gff_maker,bypass
---read_gff_maker_3UTR,gff_maker
---read_gff_maker_5UTR,gff_maker
--read_gff_maker_exon,gff_maker
</workflow>
```
##### `metadata`
The metadata part aims to provide the necessary metadata to create an EMBL entry. If plugins require specific metadata and cannot find it they will raise an exception. Metadata should be writtent as follow: `<key>:<value>`.

##### `plugins`
The plugins part aims to instantiate the various plugins and link it to a specific key. Plugins should be declared as follow: `<key>:<plugin_name>,<plugin_package>`. It is important to declare plugins required to load files as well as thos required in the workflow description.

##### `handles`
The handles part aims to load the various files with a plugin and link them with specific keys. Handles should be declared as follow: `<key>:<plugin_key>,<path_to_file>`

##### `workflow`
The workflow part aims to describe how the plugins will be called. Each task is prefixed with a certain amount of `-` which represents the level of the task. A task of level 3 will be nested in a task of level 2 which will also be nested in a task of level 1, etc. Tasks are then declared as follow: `<plugin_key>,<handle_key>,<decoration>`. At this time we have three types of decoration: `default`, `verify` and `bypass`. Each task where the decoration is not specified will be declared as `default` tasks. If a verification task fails, only the function declared as `bypass` will be called.

### Run EMBLish:
```
python3 EMBLish/main.py config.info
```

## File Validation
The output can be validated using the ENA validator found [here](https://search.maven.org/artifact/uk.ac.ebi.ena.sequence/embl-api-validator)

