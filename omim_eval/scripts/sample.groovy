def entries = new File('data/pp_conv.tsv').text.split('\n')

def r = new Random(42069)

def train = []
def test = []

entries.each {
  def g = r.nextInt(2)
  if(g == 0) {
    train << it
  } else {
    test << it
  }
}

new File('data/pp_conv_train.tsv').text = train.join('\n')
new File('data/pp_conv_test.tsv').text = test.join('\n')
