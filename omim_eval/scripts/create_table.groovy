@Grab('org.apache.commons:commons-csv:1.2')
import org.apache.commons.csv.CSVParser;
import static org.apache.commons.csv.CSVFormat.*;

def DATA_PATH = 'results'
def HEADER = '''
    \\begin{tabular}{p{8cm}|l|l|l|l|l|l}
     & \\multicolumn{3}{c|}{Binomial} & \\multicolumn{2}{c}{Fisher} \\\\
     & zscore & OR & p & zscore & OR & p \\\\
     \\hline \\hline
'''
def TAILER = '\\end{tabular}'
def dataDir = new File(DATA_PATH)

def handle = { res, test, record ->
  def iri = record['TermID']
  if(Float.parseFloat(record['adjp']) > 0.05) { return; }
  if(!res.containsKey(iri)) {
    res[iri] = [:]
  }
  res[iri][test] = record
}
def getCol = { v, test -> 
  if(v.containsKey(test)) {
    def o = v[test]
    "${o['zscore']} & ${o['or']} & ${o['adjp']}"
  } else {
    '& &' 
  }
}
def printTable = { pr ->
  def pout = HEADER + '\n'
  pr.sort { k -> 
    v = k.getValue()
    v.containsKey('binomial') ? -(Float.parseFloat(v['binomial']['zscore'])) : -(Float.parseFloat(v['fisher']['zscore']))
  }.each { k, v ->
    def infoRec = v.containsKey('binomial') ? v['binomial'] : v['fisher']
    pout += "${infoRec['name']} (\\tt ${infoRec['TermID']}) & "
    pout += getCol(v, 'binomial') + ' & '
    pout += getCol(v, 'fisher')
    pout += '\\\\ \n'
  }
  pout += TAILER
  println pout
}

// Phenopacket Results
def pr = [:]
new File('data/binomial.sigResults_enriched.tsv').withReader { reader ->
  CSVParser csv = new CSVParser(reader, TDF.withHeader())
  for(record in csv.iterator()) {
    handle(pr, 'binomial', record)
  }
}
new File('data/fisher.sigResults_enriched.tsv').withReader { reader ->
  CSVParser csv = new CSVParser(reader, TDF.withHeader())
  for(record in csv.iterator()) {
    handle(pr, 'fisher', record)
  }
}
printTable(pr)
