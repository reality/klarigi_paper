@Grab(group='com.github.sharispe', module='slib-sml', version='0.9.1')

import groovyx.gpars.*
import org.codehaus.gpars.*
import java.util.concurrent.*
 
import slib.graph.algo.utils.GAction;
import slib.graph.algo.utils.GActionType;

import org.openrdf.model.URI;
import slib.graph.algo.accessor.GraphAccessor;
//import slib.graph.algo.extraction.utils.GAction;
//import slib.graph.algo.extraction.utils.GActionType;
import slib.graph.algo.validator.dag.ValidatorDAG;
import slib.graph.io.conf.GDataConf;
import slib.graph.io.conf.GraphConf;
import slib.graph.io.loader.GraphLoaderGeneric;
import slib.graph.io.util.GFormat;
import slib.graph.model.graph.G;
import slib.graph.model.impl.graph.memory.GraphMemory;
import slib.graph.model.impl.repo.URIFactoryMemory;
import slib.graph.model.repo.URIFactory;
import slib.sml.sm.core.engine.SM_Engine;
import slib.sml.sm.core.metrics.ic.utils.IC_Conf_Topo;
import slib.sml.sm.core.metrics.ic.utils.ICconf;
import slib.graph.algo.extraction.utils.*
import slib.sglib.io.loader.*
import slib.sml.sm.core.metrics.ic.utils.*
import slib.sml.sm.core.utils.*
import slib.sglib.io.loader.bio.obo.*
import org.openrdf.model.URI
import slib.graph.algo.extraction.rvf.instances.*
import slib.sglib.algo.graph.utils.*
import slib.utils.impl.Timer
import slib.graph.algo.extraction.utils.*
import slib.graph.model.graph.*
import slib.graph.model.repo.*
import slib.graph.model.impl.graph.memory.*
import slib.sml.sm.core.engine.*
import slib.graph.io.conf.*
import slib.graph.model.impl.graph.elements.*
import slib.graph.algo.extraction.rvf.instances.impl.*
import slib.graph.model.impl.repo.*
import slib.graph.io.util.*
import slib.graph.io.loader.*

import slib.sml.sm.core.metrics.ic.utils.IC_Conf_Corpus;
import slib.sml.sm.core.utils.SMConstants;
import slib.sml.sm.core.utils.SMconf;
import slib.utils.ex.SLIB_Exception;
import slib.utils.impl.Timer;

//def omim = new File('hpo_mimc_anns.txt').text.split('\n').collect { it.split('\t')[0] }
def multivariate = [:]
new File('data/associations/multivariate_klarigi.tsv').splitEachLine('\t') {
  if(!multivariate.containsKey(it[1])) { multivariate[it[1]] = [] }
  multivariate[it[1]] << it[0]
}

def univariate = [:]
new File('data/associations/univariate_klarigi.tsv').splitEachLine('\t') {
  if(!univariate.containsKey(it[1])) { univariate[it[1]] = [] }
  univariate[it[1]] << it[0]
}
/*def univariate = [:]
new File('data/associations/univariate_klarigi.tsv').splitEachLine {
  if(!multivariate.containsKey(it[1])) { multivariate[it[1]] = [] }
  multivariate[it[1]] << it[0]
}*/

def enrichment = [:]
new File('data/associations/enrichment.tsv').splitEachLine('\t') {
  if(!enrichment.containsKey(it[1])) { enrichment[it[1]] = [] }
  enrichment[it[1]] << it[0]
}

println 'writing the annotation file now'
def sWriter = new BufferedWriter(new FileWriter('data/annot.tsv'))
def oo = ''
def y = 0


/*dMap.each { a, b ->
  println "(${++y}/${aMap.size()})"
  if(!b.any{ it.indexOf('txt') != -1}) {
    sWriter.write('http://reality.rehab/ptvis/' + a + '\t' + b.join(';') + '\n')
  }
}
sWriter.flush()
sWriter.close()
println 'done'*/

URIFactory factory = URIFactoryMemory.getSingleton()

def ontoFile = 'data/hp.owl'
def graphURI = factory.getURI('http://HP/')
factory.loadNamespacePrefix("HP", graphURI.toString());

G graph = new GraphMemory(graphURI)

def dataConf = new GDataConf(GFormat.RDF_XML, ontoFile)
def actionRerootConf = new GAction(GActionType.REROOTING)
actionRerootConf.addParameter("root_uri", "http://purl.obolibrary.org/obo/HP_0000001"); // phenotypic abnormality
//actionRerootConf.addParameter("root_uri", "DOID:4"); // phenotypic abnormality

def gConf = new GraphConf()
gConf.addGDataConf(dataConf)
gConf.addGAction(actionRerootConf)
//def gConf = new GraphConf()
//gConf.addGDataConf(dataConf)
def annot = 'data/train.tsv'
gConf.addGDataConf(new GDataConf(GFormat.TSV_ANNOT, annot));

GraphLoaderGeneric.load(gConf, graph)

def roots = new ValidatorDAG().getTaxonomicRoots(graph)

def icConf = new IC_Conf_Corpus(SMConstants.FLAG_IC_ANNOT_RESNIK_1995)
//def icConf = new IC_Conf_Topo(SMConstants.FLAG_ICI_ZHOU_2008)
//def icConf = new IC_Conf_Topo(SMConstants.FLAG_ICI_SANCHEZ_2011)
def smConfPairwise = new SMconf(SMConstants.FLAG_SIM_PAIRWISE_DAG_NODE_RESNIK_1995, icConf)
//def smConfPairwise = new SMconf(SMConstants.FLAG_SIM_PAIRWISE_DAG_NODE_LIN_1998, icConf)

def smConfGroupwise = new SMconf(SMConstants.FLAG_SIM_GROUPWISE_AVERAGE, icConf)
// uncomment to do the bma
//def smConfGroupwise = new SMconf(SMConstants.FLAG_SIM_GROUPWISE_BMA, icConf)

//def smConfGroupwise = new SMconf(SMConstants.FLAG_SIM_GROUPWISE_AVERAGE, icConf)
// FLAG_SIM_GROUPWISE_AVERAGE_NORMALIZED_GOSIM

//def smConfPairwise = new SMconf(SMConstants.FLAG_SIM_PAIRWISE_DAG_NODE_JIANG_CONRATH_1997_NORM , icConf)


def z = 0

def engine = new SM_Engine(graph)

def getURIfromTerm = { term ->
  term = term.tokenize(':')
  return factory.getURI("http://purl.obolibrary.org/obo/HP_" + term[1])
}

def sims = []
def out = ["pesim\tpnsim\tmatch1\tmatch2"]
new File('data/test.tsv').splitEachLine('\t') {
def pesim = engine.compare(smConfGroupwise, smConfPairwise, enrichment['pulmonary embolism'].collect { getURIfromTerm(it) }.toSet(), it[1].split(';').collect { getURIfromTerm(it) }.findAll { graph.containsVertex(it) }.toSet())
def pnsim = engine.compare(smConfGroupwise, smConfPairwise, enrichment['pneumonia'].collect { getURIfromTerm(it) }.toSet(), it[1].split(';').collect { getURIfromTerm(it) }.findAll { graph.containsVertex(it) }.toSet())
def match1 = it[2] == 'pulmonary embolism'
def match2 = it[2] == 'pneumonia'
  out << "$pesim\t$pnsim\t$match1\t$match2"
}
new File('data/sim/enrichment_sim_class.tsv').text = out.join('\n')

out = ["pesim\tpnsim\tmatch1\tmatch2"]

new File('data/test.tsv').splitEachLine('\t') {
def pesim = engine.compare(smConfGroupwise, smConfPairwise, multivariate['pulmonary embolism'].collect { getURIfromTerm(it) }.toSet(), it[1].split(';').collect { getURIfromTerm(it) }.findAll { graph.containsVertex(it) }.toSet())
def pnsim = engine.compare(smConfGroupwise, smConfPairwise, multivariate['pneumonia'].collect { getURIfromTerm(it) }.toSet(), it[1].split(';').collect { getURIfromTerm(it) }.findAll { graph.containsVertex(it) }.toSet())
def match1 = it[2] == 'pulmonary embolism'
def match2 = it[2] == 'pneumonia'
  out << "$pesim\t$pnsim\t$match1\t$match2"
}
new File('data/sim/multivariate_sim_class.tsv').text = out.join('\n')

out = ["pesim\tpnsim\tmatch1\tmatch2"]

new File('data/test.tsv').splitEachLine('\t') {
def pesim = engine.compare(smConfGroupwise, smConfPairwise, univariate['pulmonary embolism'].collect { getURIfromTerm(it) }.toSet(), it[1].split(';').collect { getURIfromTerm(it) }.findAll { graph.containsVertex(it) }.toSet())
def pnsim = engine.compare(smConfGroupwise, smConfPairwise, univariate['pneumonia'].collect { getURIfromTerm(it) }.toSet(), it[1].split(';').collect { getURIfromTerm(it) }.findAll { graph.containsVertex(it) }.toSet())
def match1 = it[2] == 'pulmonary embolism'
def match2 = it[2] == 'pneumonia'
  out << "$pesim\t$pnsim\t$match1\t$match2"
}
new File('data/sim/univariate_sim_class.tsv').text = out.join('\n')


println 'done'
