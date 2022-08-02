def entries = new File('./data/pulm_pneum_primaries_nodubs.tsv').text.split('\n')

def r = new Random(42069)

def train = []
def test = []

entries.each {
  def g = r.nextInt(5)
  if(g == 0) {
    test << it
  } else {
    train << it
  }
}

new File('data/train.tsv').text = train.join('\n')
new File('data/test.tsv').text = test.join('\n')
