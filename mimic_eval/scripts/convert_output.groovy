// this is so stupid
def tables = new File('klarigi_outputs/overall/main_ecm_withp.out').text.split('tabular')

new File('data/associations/univariate_pn_klarigi.tsv').text = (tables[1] =~ /.*(HP\:.......).*/).collect { it[1] + "\tpneumonia" }.join('\n')
new File('data/associations/univariate_pe_klarigi.tsv').text = (tables[3] =~ /.*(HP\:.......).*/).collect { it[1] + "\tpulmonary embolism" }.join('\n')

new File('data/associations/multivariate_pn_klarigi.tsv').text = (tables[5] =~ /.*(HP\:.......).*/).collect { it[1] + "\tpneumonia" }.join('\n')
new File('data/associations/multivariate_pe_klarigi.tsv').text = (tables[7] =~ /.*(HP\:.......).*/).collect { it[1] + "\tpulmonary embolism" }.join('\n')
