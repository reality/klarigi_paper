time klarigi --data data/train.tsv -o data/hp.owl --verbose --threads 6 --output-type=latex --output-scores -ecm --reclassify --resnik-ic --min-power=0.1 --min-exclusion=0.125 --perms=2000 > klarigi_outputs/main_ecm_withp.out
#time klarigi --data data/train.tsv -o ../hp.owl --verbose --threads 6 --output-type=latex --output-scores -ecm --reclassify --resnik-ic --min-power=0.1 --min-exclusion=0.125 > klarigi_outputs/main_ecm_nop.out
rm data/scores_examples/*
mv HP_*.txt data/scores_examples
groovy scripts/convert_output.groovy
