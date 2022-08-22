@Grab('org.apache.commons:commons-csv:1.2')
import org.apache.commons.csv.CSVParser;
import static org.apache.commons.csv.CSVFormat.*;

def out = []

new File('enrichment/binomial.sigResults_enrichedPulmonaryEmbolism.tsv').withReader { reader ->
  CSVParser csv = new CSVParser(reader, TDF.withHeader())
  for(record in csv.iterator()) {
    out << record['TermID'] + '\tpulmonary embolism'
  }
}

new File('enrichment/fisher.sigResults_enrichedPulmonaryEmbolism.tsv').withReader { reader ->
  CSVParser csv = new CSVParser(reader, TDF.withHeader())
  for(record in csv.iterator()) {
    out << record['TermID'] + '\tpulmonary embolism'
  }
}

new File('enrichment/binomial.sigResults_enrichedPneumonia.tsv').withReader { reader ->
  CSVParser csv = new CSVParser(reader, TDF.withHeader())
  for(record in csv.iterator()) {
    out << record['TermID'] + '\tpneumonia'
  }
}
new File('enrichment/fisher.sigResults_enrichedPneumonia.tsv').withReader { reader ->
  CSVParser csv = new CSVParser(reader, TDF.withHeader())
  for(record in csv.iterator()) {
    out << record['TermID'] + '\tpneumonia'
  }
}

out = out.unique(false)

new File('data/associations/enrichment.tsv').text = out.join('\n')
