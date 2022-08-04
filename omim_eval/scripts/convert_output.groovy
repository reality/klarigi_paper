// this is so stupid
def tables = new File('klarigi_outputs/ecm_withp.out').text.split('tabular')

def uv = (tables[1] =~ /.*(HP\:.......).*/).collect { it[1] + "\tOMIM:616900" }.join('\n')
new File('data/associations/univariate_klarigi.tsv').text = uv

def mv = (tables[3] =~ /.*(HP\:.......).*/).collect { it[1] + "\tOMIM:616900" }.join('\n')
new File('data/associations/multivariate_klarigi.tsv').text = mv

