which klarigi
time klarigi --data data/train.tsv -o data/hp.owl --verbose --threads 6 --resnik-ic --output-type=latex --classify data/test.tsv --classify-with-variables=data/associations/multivariate_klarigi.tsv --output-classification-scores > klarigi_outputs/test_klarigi_multivariate.out
#time klarigi --data data/train.tsv --classify data/test.tsv -o data/hp.owl --verbose --threads 6 --output-type=latex --reclassify --classify-with-variables=data/associations/univariate_klarigi.tsv --output-classification-scores > klarigi_outputs/test_klarigi_univariate.out
