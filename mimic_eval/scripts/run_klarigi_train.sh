which klarigi
time klarigi --data data/train.tsv -o data/hp.owl --verbose --threads 6 --output-type=latex --output-scores --reclassify --resnik-ic --min-inclusion=0.125 > klarigi_outputs/main_ecm_nop.out
#time klarigi --data data/train.tsv -o data/hp.owl --verbose --threads 6 --output-type=latex --output-scores --reclassify --resnik-ic --perms=2000 --min-inclusion=0.125 > klarigi_outputs/main_ecm_withp.out
#rm data/scores_examples/*
#mv HP_*.txt data/scores_examples
groovy scripts/convert_output.groovy
