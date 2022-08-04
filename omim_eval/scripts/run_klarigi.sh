klarigi --data data/phenopackets -pp -o data/hp.owl --verbose --threads 6 --output-scores --group OMIM:616900 --output-type=latex --resnik-ic --min-power=0.2 --reclassify -ecm --perms=2000 > klarigi_outputs/ecm_withp.out 
klarigi --data data/phenopackets -pp -o data/hp.owl --verbose --threads 6 --group OMIM:616900 --output-type=latex --resnik-ic --min-power=0.2 --reclassify -ecm > klarigi_outputs/ecm_nop.out 
rm data/scores/*
mv HP_*.txt data/scores
