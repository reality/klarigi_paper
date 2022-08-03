// this is so stupid
def tables = new File('klarigi_outputs/overall/main_ecm_withp.out').text.split('tabular')

def uv_pe = (tables[1] =~ /.*(HP\:.......).*/).collect { it[1] + "\tpulmonary embolism" }.join('\n')
def uv_pn = (tables[3] =~ /.*(HP\:.......).*/).collect { it[1] + "\tpneumonia" }.join('\n')
new File('data/associations/univariate_pn_klarigi.tsv').text = uv_pn
new File('data/associations/univariate_pe_klarigi.tsv').text = uv_pe
new File('data/associations/univariate_klarigi.tsv').text = uv_pe + '\n' + uv_pn

def mv_pn = (tables[5] =~ /.*(HP\:.......).*/).collect { it[1] + "\tpneumonia" }.join('\n')
def mv_pe = (tables[7] =~ /.*(HP\:.......).*/).collect { it[1] + "\tpulmonary embolism" }.join('\n')
new File('data/associations/multivariate_pn_klarigi.tsv').text = mv_pn
new File('data/associations/multivariate_pe_klarigi.tsv').text = mv_pe
new File('data/associations/multivariate_klarigi.tsv').text = mv_pe + '\n' + mv_pn

