def entries = new File('./data/pulm_pneum_primaries_nodubs.tsv').text.split('\n')

def STRIP = [ "HP:0100598", "HP:0002090", "HP:0011949", "HP:0011951", "HP:0002100", "HP:0011952" ]

def r = new Random(42069)

def train = []
def test = []

entries.each {
  def g = r.nextInt(5)

  it = it.split('\t')
  it[1] = it[1].split(';').collect { 'HP:' + it.split('_').last() }.findAll { i -> !STRIP.contains(i) }.join(';')
  it = it.join('\t')

  if(g == 0) {
    test << it
  } else {
    train << it
  }
}

new File('data/train.tsv').text = train.join('\n')
new File('data/test.tsv').text = test.join('\n')
