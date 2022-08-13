which klarigi
klarigi --data data/phenopackets -pp -o data/hp.owl --output-scores --verbose --threads 6 --group OMIM:616900 --output-type=latex --resnik-ic --min-r-score=0.2 --reclassify > klarigi_outputs/ecm_nop.out 
#klarigi --data data/phenopackets -pp -o data/hp.owl --verbose --threads 6 --output-scores --group OMIM:616900 --output-type=latex --resnik-ic --min-r-score=0.2 --reclassify --perms=2000 > klarigi_outputs/ecm_withp.out 
#rm data/scores/*
#mv HP_*.txt data/scores
