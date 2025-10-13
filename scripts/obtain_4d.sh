bin/phast-1.3/bin/msa_view  splitmaf/simHuman1chr7.maf --4d --features simchr7.gff > 4d-codons.ss
bin/phast-1.3/bin/msa_view 4d-codons.ss --in-format SS  --tuple-size 1 > 4d-sites.fas
