@Grapes([
    @Grab(group='net.sourceforge.owlapi', module='owlapi-api', version='5.1.20'),
    @Grab(group='net.sourceforge.owlapi', module='owlapi-apibinding', version='5.1.20'),
    @Grab(group='net.sourceforge.owlapi', module='owlapi-impl', version='5.1.20'),
    @Grab(group='net.sourceforge.owlapi', module='owlapi-parsers', version='5.1.20'),
    @Grab(group='net.sourceforge.owlapi', module='owlapi-distribution', version='5.1.20'),

  @GrabResolver(name='sonatype-nexus-snapshots', root='https://oss.sonatype.org/service/local/repositories/snapshots/content/'),
   // @Grab('org.semanticweb.elk:elk-hpReasoner:0.5.0-SNAPSHOT'),
   // @Grab('org.semanticweb.elk:elk-owl-implementation:0.5.0-SNAPSHOT'),
    @Grab('au.csiro:elk-owlapi5:0.5.0'),
  
    @GrabConfig(systemClassLoader=true)
])

import org.semanticweb.owlapi.model.IRI
import org.semanticweb.owlapi.model.parameters.*
import org.semanticweb.elk.owlapi.*
import org.semanticweb.elk.reasoner.config.*
import org.semanticweb.owlapi.apibinding.OWLManager
import org.semanticweb.owlapi.reasoner.*
import org.semanticweb.owlapi.vocab.OWLRDFVocabulary
import org.semanticweb.owlapi.model.*
import org.semanticweb.owlapi.io.*
import org.semanticweb.owlapi.owllink.*
import org.semanticweb.owlapi.util.*
import org.semanticweb.owlapi.search.*
import org.semanticweb.owlapi.manchestersyntax.renderer.*
import org.semanticweb.owlapi.reasoner.structural.*
import org.semanticweb.elk.reasoner.config.*
import org.semanticweb.owlapi.apibinding.*
import org.semanticweb.owlapi.reasoner.*
import java.util.concurrent.*
import java.util.concurrent.atomic.*
import groovyx.gpars.*
import org.codehaus.gpars.*

def entries = new File('./data/pulm_pneum_primaries_nodubs.tsv').text.split('\n')

def STRIP = [ "HP:0002090", "HP:0002204", "HP:0001977" ]

def manager = OWLManager.createOWLOntologyManager()
def fac = manager.getOWLDataFactory()
def config = new SimpleConfiguration()
def elkFactory = new ElkReasonerFactory() // cute

def hpOntology = manager.loadOntologyFromOntologyDocument(new File("data/hp.owl"))
def hpReasoner = elkFactory.createReasoner(hpOntology, config)

def allStrip = STRIP.collect { hp ->
  def ce = fac.getOWLClass(IRI.create('http://purl.obolibrary.org/obo/' + hp.replace(':', '_')))
  hpReasoner.getSubClasses(ce, false).collect { it.getRepresentativeElement().getIRI().toString() }.collect { it.split('/').last().replace('_', ':') }.flatten().unique(false).findAll { it != 'owl#Nothing' } + hp
 }.flatten()

println "Full list of strips: ${allStrip.join(',')}"

def r = new Random(1242069420)

def train = []
def test = []

entries.each {
  def g = r.nextInt(5)

  it = it.split('\t')
  def og = it[1].split(';').collect { 'HP:' + it.split('_').last() }
  it[1] = og.findAll { i -> !allStrip.contains(i) }.join(';')
  it = it.join('\t')

  if(g == 0) {
    test << it
  } else {
    train << it
  }
}

new File('data/train.tsv').text = train.join('\n')
new File('data/test.tsv').text = test.join('\n')
