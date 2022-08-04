def currentInfo = []
def out = [ 
  [ "term_id", "term_name", "is_a", "is_a2", "is_a3", "is_a4" ]
]
new File('res.txt').splitEachLine(':') {
  it[1] = it.subList(1, it.size()).join(':').trim()
  it[1] = it[1].replaceAll(',','')
  if(it[0] == 'id') {
    while(currentInfo.size() < 6) {
      currentInfo << ""
    }

    out << currentInfo
    currentInfo = []
  }

  currentInfo << it[1]
}

out.each {
  println it.join(',')
}
