time klarigi --data data/train.tsv -o ../hp.owl --verbose --threads 6 --output-type=latex --output-scores -ecm --reclassify --min-power=0.125 --resnik-ic --perms=2000 > klarigi_outputs/overall/main_ecm_withp.out
mv HP_*.txt data/scores_examples
