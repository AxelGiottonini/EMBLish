# EMBLish

## Config file example
```
python3 main.py config.info
```

```
cd
|----config.info
|----sequences.fasta
|----data.gff
|----anno.out
```


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
---verify_gff_maker_CDS,NF,verify
---read_gff_maker_misc_feature,gff_maker,bypass
---read_gff_maker_3UTR,gff_maker
---read_gff_maker_5UTR,gff_maker
--read_gff_maker_exon,gff_maker
</workflow>

```
