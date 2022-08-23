time klarigi --data data/train.tsv -o data/hp.owl --verbose --threads 6 --output-type=latex --reclassify --resnik-ic --group "pulmonary embolism" -egl --min-inclusion=0.3 > klarigi_outputs/pe_egl.out
time klarigi --data data/train.tsv -o data/hp.owl --verbose --threads 6 --output-type=latex --reclassify --resnik-ic --group pneumonia -egl --max-inclusion=0.7 > klarigi_outputs/pn_egl.out
