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
/*def pr = [:]
new File(dataDir, 'enrichedBinomialBonferroniPhenoPacket.tsv').withReader { reader ->
  CSVParser csv = new CSVParser(reader, TDF.withHeader())
  for(record in csv.iterator()) {
    handle(pr, 'binomial', record)
  }
}
new File(dataDir, 'enrichedFisherBonferroniPhenoPacket.tsv').withReader { reader ->
  CSVParser csv = new CSVParser(reader, TDF.withHeader())
  for(record in csv.iterator()) {
    handle(pr, 'fisher', record)
  }
}
printTable(pr)*/

//new File('omimclas.tsv').text = pr.collect { k, v -> "$k\tOMIM:616900" }.join('\n')

def PE_PN_REMOVE = ['HP:0100598', 'HP:0100598', 'HP:0002090', 'HP:0011949', 'HP:0011951', 'HP:0011952', 'HP:0002100', 'HP:0006532', 'HP:0033214']

def pe = [:]
new File('data/enrichment/binomial.sigResults_enrichedPulmonaryEmbolism.tsv').withReader { reader ->
  CSVParser csv = new CSVParser(reader, TDF.withHeader())
  for(record in csv.iterator()) {
    if(!PE_PN_REMOVE.contains(record['TermID'])) {
      handle(pe, 'binomial', record)
    } else { println ':)'}
  }
}
new File('data/enrichment/fisher.sigResults_enrichedPulmonaryEmbolism.tsv').withReader { reader ->
  CSVParser csv = new CSVParser(reader, TDF.withHeader())
  for(record in csv.iterator()) {
    if(!PE_PN_REMOVE.contains(record['TermID'])) {
      handle(pe, 'fisher', record)
    } else { println ':)'}
  }
}
printTable(pe)

def pn = [:]
new File('data/enrichment/binomial.sigResults_enrichedPneumonia.tsv').withReader { reader ->
  CSVParser csv = new CSVParser(reader, TDF.withHeader())
  for(record in csv.iterator()) {
    if(!PE_PN_REMOVE.contains(record['TermID'])) {
      handle(pn, 'binomial', record)
    }
  }
}
new File('data/enrichment/fisher.sigResults_enrichedPneumonia.tsv').withReader { reader ->
  CSVParser csv = new CSVParser(reader, TDF.withHeader())
  for(record in csv.iterator()) {
    if(!PE_PN_REMOVE.contains(record['TermID'])) {
      handle(pn, 'fisher', record)
    }
  }
}
printTable(pn)

