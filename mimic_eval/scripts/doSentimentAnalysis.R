# Sentiment analysis for klarigi paper
# @ John Williams - j.a.williams@bham.ac.uk ; johnwilliams@post.harvard.edu
# main functions credit XGR, igraph folks 
install.packages("BiocManager")
BiocManager::install("gghex", dependencies=T)
BiocManager::install("hfang-bristol/XGR", dependencies=T)
## reload the installed package
detach(package:XGR, unload=T)

# check packages and install if missing
using<-function(...) {
  libs<-unlist(list(...))
  req<-unlist(lapply(libs,require,character.only=TRUE))
  need<-libs[req==FALSE]
  if(length(need)>0){ 
    install.packages(need)
    lapply(need,require,character.only=TRUE)
  }
}
using("tidyverse", "data.table","igraph","qdapTools","Matrix")

# note - can get away without tidyverse but individual packages like
library(XGR)
# dplyr and magrittr would need to be installed and John is lazy so ... 

### also the function above should load all libraies needed but
# if a function not found just load one of them manually 

# working directory contains hp.obo, pulm_pneum_primaries_nodubs.csv, pp_conv.csv 

setwd("enrichment/")

# John was lazy so quick munge of obo relations to res.txt
system("perl ./extract_rels.perl > res.txt")
system("groovy transform.groovy > mat.lst")

############ perl extract_rels looks like:
#	
#	#!/usr/bin/perl
#	use strict; use warnings; use autodie;

#	open my $handle, '<', qq(hp.obo);
#	chomp(my @lines = <$handle>);
#	close $handle;

#	foreach (@lines) {

#	    if ($_=~ /^id/) {
#	        print(qq($_\n));
#	    }
#	        if ($_=~ /^name/) {
#	        print(qq($_\n));
#	    }
#	        if ($_=~ /^is_a/) {
#	        $_ =~ s/ !.*//g;
#	        print(qq($_\n));
#	    }
#	}




# pick your poison --- load either experiment
#dat <- data.table::fread("pulm_pneum_primaries_nodubs.csv", data.table = F)
#dat <- data.table::fread("pp_conv.tsv")

dat <- data.table::fread("../data/train.tsv", data.table = F)
# for XGR to work well (which can handle DAG-based p-val correction) 
# we need to check structure of thier pre-loaded HPO annotations 
### note - probably subset theirs to HPPA and not having the pseudo-root nodes
### is what caused the code to want to fail on phenopackets 

gHPO <- xRDataLoader("ig.HPPA")
gHPO %>% summary()

V(gHPO)$name
V(gHPO)$term_id
V(gHPO)$term_name
V(gHPO)$term_distance
E(gHPO)$relation %>% table()
is.directed(gHPO)
gHPO

# now remove gHPO cuz bleh 
rm(gHPO)

# red in matrix, remove obsoltete classes 


#### now would be where to remove list of non- HPPA classes but John didn't do that yet

# read
#myMat <- readLines("mat.lst")
# rm obsolete
#myMat <- myMat[!grepl("obsolete",myMat)]
# split parent classes by ;
#myMat <- str_split(myMat, ";")
#myMat <- plyr::ldply(myMat, rbind)
# check dimensions and head  to know how many columns to add 
#dim(myMat)
#head(myMat)
# remove first row (pseudo column names)
#myMat <- myMat[-1,]

# add column names for each is_a relation (parent/child)

myMat <- read.csv("mat.lst")
colnames(myMat) <- c("term_id","term_name","is_a", "is_a2", "is_a3","is_a4")
myMat <- as.matrix(myMat)

# create HP dag in list form of parent, child
L1 <- myMat[,c(3,1)]
L2 <- myMat[,c(4,1)]
L3 <- myMat[,c(5,1)]
L4 <- myMat[,c(6,1)]
colnames(L1) <- c("Parent","Child")
colnames(L2) <- c("Parent","Child")
colnames(L3) <- c("Parent","Child")
colnames(L4) <- c("Parent","Child")
# bind all together into L1 and remove others
L1 <- rbind(L1,L2)
L1 <- rbind(L1,L3)
L1 <- rbind(L1,L4)
rm(L2,L3,L4)
# remove missing parent info 
L1 <- L1[complete.cases(L1),]

L1 <- L1[L1[,1] != "",]

# create node list using first two columns of matrix as start 
myNodeList <- myMat[,1:2]
myNodeList <- as.data.frame(myMat)
# create HPO graph
HPGraph <- graph_from_edgelist(as.matrix(L1), directed = TRUE) %>% igraph::simplify()
# check that the graph worked
HPGraph
# subset to remove any terms excluded earlier 
myNodeList <- myNodeList[myNodeList$term_id %in% names(V(HPGraph)),]
# now the fun part -- add vertex and edge attributes, hopefully self explanitory 

HPGraph <- set_vertex_attr(HPGraph, "term_id", index = myNodeList$term_id, myNodeList$term_id)
HPGraph <- set_vertex_attr(HPGraph, "term_name", index = myNodeList$term_id, myNodeList$term_name)

HPGraph <- set_edge_attr(HPGraph, "relation",value = "is_a")

# now we need to add path distance for DAG stuff --- so use shortest_path 
# from igraph library to build up 
mypaths <- shortest_paths(HPGraph, from = 1, mode = "out", output = "both")
# helper vars all related to getting out distances and names ...  
myLengths <- mypaths[[1]]
myLengthsName <- lapply(myLengths, FUN = function(x) {tail(x, n = 1)}) %>% unlist() %>% names()
myLenghtsLength <- lapply(myLengths, length) %>% unlist()
myLenghts <- myLenghtsLength
names(myLengths) <- myLengthsName
# offset by one to agree with the pre-computed gHPO that XGR uses (so path from self to self is 0, not 1)
myLenghts <- myLenghts - 1
names(myLenghts) <- myLengthsName
# make a dataframe out fo this 
LDf <- data.frame(term_id = names(myLengths), term_distance = myLenghts)
# get term IDs 'in order' so can sort lenths correctly 
termDistances <- data.frame(V(HPGraph)$term_id)
# data munging ... 
termDistances <- dplyr::left_join(termDistances,LDf, by = c("V.HPGraph..term_id" = "term_id"))
# set unconnected nodes to distance 0 
termDistances$term_distance[termDistances$term_distance == -1] <- 0
# make sure no extra distances added cuz of internal checks in XGR 
termDistances <- termDistances[termDistances$V.HPGraph..term_id %in% V(HPGraph)$term_id,]
termDistances <- unique(termDistances)
# final dummy thing then join all up to add the vertex attribute of term distance 
myNodeList <- dplyr::left_join(myNodeList, termDistances, by = c("term_id" = "V.HPGraph..term_id"))
HPGraph <- set_vertex_attr(HPGraph, "term_distance", index = myNodeList$term_id, myNodeList$term_distance) 


##### so here, start with dat (the data being tested) and make groups
#### the group name 'pneumoniaGroup' left over from doing the non-phenopackets example
#### but phenopacets used here I think :) 
colnames(dat) <- c("Transaction","Terms","Group")
pneumoniaGroup <- dat$Transaction[dat$Group == "pneumonia"]
pulmonaryEmbolismGroup <- dat$Transaction[dat$Group == "pulmonary embolism"]

# gsub IRI 
dat$Terms <- gsub("http://purl.obolibrary.org/obo/","",dat$Terms)
dat$Terms <- gsub("_",":",dat$Terms)
datList <- strsplit(dat$Terms,";")
names(datList) <- dat$Transaction

# change format of data for graph work then coerese to matrix 
d1 <- qdapTools::mtabulate(setNames(strsplit(as.character(dat$Terms), ";", 
                                  fixed=T), dat$Transaction))
d1 <- Matrix::Matrix(as.matrix(d1),sparse = T)



##### example on how to test ... again not duplicated for each case cuz John lazy ... 
##### what will change is the data and the test (to fisher)
##### here we do bonferroni corretion and correct pvals conservatively for DAG structure
##### so we can say we're using state of art such and such ... 
##### testing background == every annotated patient 

##### go to line 215-ish and load xEnricher1 function ... 

resPneumonia <- xEnricher1(data = pneumoniaGroup, # change this when switching things up
                          annotation = d1, # and this make new for each pheno dataset 
                          g = HPGraph,
                          min.overlap = 3,
                          test = "binomial", # change to fisher as well 
                          p.tail = "one-tail",
                          background = rownames(d1),
                          p.adjust.method = "bonferroni",
                          verbose = TRUE,
                          true.path.rule = TRUE,
                          ontology.algorithm = "pc")

#### coerse results into dataframe and add the tested group in case 
#### joined up with other resuls 
resPneumoniaDF <- xEnrichViewer(resPneumonia, top_num = length(resPneumonia$adjp))
resPneumoniaDF$groupTested <- "pneumonia"

#### example of changing it up assuming switched d1 
respulmonaryEmbolismGroup <- xEnricher1(data = pulmonaryEmbolismGroup, 
                           annotation = d1,
                           g = HPGraph,
                           min.overlap = 3,
                           test = "binomial",
                           p.tail = "one-tail",
                           background = rownames(d1),
                           p.adjust.method = "bonferroni",
                           verbose = TRUE,
                           true.path.rule = TRUE,
                           ontology.algorithm = "pc")


respulmonaryEmbolismGroupDF <- xEnrichViewer(respulmonaryEmbolismGroup, top_num = length(respulmonaryEmbolismGroup$adjp))
respulmonaryEmbolismGroupDF$groupTested <- "pulmonaryEmbolism"
### to make things better we can add the term ID from row names cuz they're stored there 
respulmonaryEmbolismGroupDF$TermID <- rownames(respulmonaryEmbolismGroupDF)
resPneumoniaDF$TermID <- rownames(resPneumoniaDF)
## and we can write things out ... 
#write.table(respulmonaryEmbolismGroupDF,"enrichedPulmonaryEmbolism.tsv",row.names = F, col.names = T, quote = F, sep = "\t")
#write.table(resPneumoniaDF,"enrichedPneumonia.tsv",row.names = F, col.names = T, quote = F, sep = "\t")

## or we can write just the adjusted corrected data 

sigResp <- respulmonaryEmbolismGroupDF[respulmonaryEmbolismGroupDF$adjp < 0.05,]
sigPneu <- resPneumoniaDF[resPneumoniaDF$adjp < 0.05,]
write.table(sigResp,"binomial.sigResults_enrichedPulmonaryEmbolism.tsv",row.names = F, col.names = T, quote = F, sep = "\t")
write.table(sigPneu,"binomial.sigResults_enrichedPneumonia.tsv",row.names = F, col.names = T, quote = F, sep = "\t")

#write.table(omim_sigPneu,"binomial.omim.dianosis.tsv",row.names = F, col.names = T, quote = F, sep = "\t")


##### below function must beloaded first ... 



respulmonaryEmbolismGroup <- xEnricher1(data = pulmonaryEmbolismGroup, 
                                        annotation = d1,
                                        g = gHPO,
                                        min.overlap = 3,
                                        test = "binomial",
                                        p.tail = "one-tail",
                                        background = rownames(d1),
                                        p.adjust.method = "BH",
                                        verbose = TRUE,
                                        true.path.rule = TRUE,
                                        ontology.algorithm = "pc")
respulmonaryEmbolismGroupDF <- xEnrichViewer(respulmonaryEmbolismGroup, top_num = length(respulmonaryEmbolismGroup$adjp))
respulmonaryEmbolismGroupDF$groupTested <- "pulmonaryEmbolism"
respulmonaryEmbolismGroupDF$TermID <- rownames(respulmonaryEmbolismGroupDF)
resPneumoniaDF$TermID <- rownames(resPneumoniaDF)


### now do phenopackets

phenopacks <- read_tsv("pp_conv (2).csv")
table(phenopacks$`OMIM:191900`) %>% sort()



















#### function modified from XGR package by John 
xEnricher1 <- function (data, annotation, g, background = NULL, size.range = c(10, 
                                                                               2000), min.overlap = 5, which.distance = NULL, test = c("fisher", 
                                                                                                                                       "hypergeo", "binomial"), background.annotatable.only = NULL, 
                        p.tail = c("one-tail", "two-tails"), p.adjust.method = c("BH", 
                                                                                 "BY", "bonferroni", "holm", "hochberg", "hommel"), ontology.algorithm = c("none", 
                                                                                                                                                           "pc", "elim", "lea"), elim.pvalue = 0.01, lea.depth = 2, 
                        path.mode = c("all_paths", "shortest_paths", "all_shortest_paths"), 
                        true.path.rule = TRUE, verbose = TRUE) 
{
  test <- match.arg(test)
  p.tail <- match.arg(p.tail)
  p.adjust.method <- match.arg(p.adjust.method)
  ontology.algorithm <- match.arg(ontology.algorithm)
  path.mode <- match.arg(path.mode)
  p.tail <- match.arg(p.tail)
  if (is.vector(data)) {
    data <- unique(data)
  }
  else {
    warnings("The input data must be a vector.\n")
    return(NULL)
  }
  if (is(annotation, "GS")) {
    originAnnos <- annotation$gs
  }
  else if (is(annotation, "list")) {
    originAnnos <- annotation
  }
  else if (is(annotation, "dgCMatrix")) {
    D <- annotation
    originAnnos <- sapply(1:ncol(D), function(j) {
      names(which(D[, j] != 0))
    })
    names(originAnnos) <- colnames(annotation)
  }
  else {
    warnings("The input annotation must be either 'GS' or 'list' or 'dgCMatrix' object.\n")
    return(NULL)
  }
  annotation <- originAnnos
  ig <- g
  if (!is(ig, "igraph")) {
    warnings("The function must apply to the 'igraph' object.\n")
    return(NULL)
  }
  else {
    if (verbose) {
      now <- Sys.time()
      message(sprintf("First, generate a subgraph induced (via '%s' mode) by the annotation data (%s) ...", 
                      path.mode, as.character(now)), appendLF = TRUE)
    }
    subg <- xDAGanno(g = ig, annotation = annotation, path.mode = path.mode, 
                     true.path.rule = true.path.rule, verbose = verbose)
    gs <- V(subg)$anno
    names(gs) <- V(subg)$name
    gs.distance <- V(subg)$term_distance
    names(gs.distance) <- V(subg)$name
  }
  if (1) {
    if (is.vector(background)) {
      background <- base::unique(background)
      background <- background[!is.null(background)]
      background <- background[!is.na(background)]
    }
    if (length(background) > 0) {
      if (1) {
        background <- base::union(background, data)
      }
      gs <- lapply(gs, function(x) {
        ind <- match(x, background)
        x[!is.na(ind)]
      })
    }
  }
  if (!is.null(which.distance) & sum(is.na(gs.distance)) == 
      0) {
    distance_filtered <- lapply(which.distance, function(x) {
      names(gs)[(gs.distance == as.integer(x))]
    })
    distance_filtered <- unlist(distance_filtered)
  }
  else {
    distance_filtered <- names(gs)
  }
  ind.distance <- match(distance_filtered, names(gs))
  gs.length <- sapply(gs, length)
  ind.length <- which(gs.length >= size.range[1] & gs.length <= 
                        size.range[2])
  ind <- intersect(ind.distance, ind.length)
  gs <- gs[ind]
  if (length(gs) == 0) {
    warnings("There are no terms being used.\n")
    return(NULL)
  }
  doFisherTest <- function(genes.group, genes.term, genes.universe) {
    genes.hit <- intersect(genes.group, genes.term)
    X <- length(genes.hit)
    K <- length(genes.group)
    M <- length(genes.term)
    N <- length(genes.universe)
    cTab <- matrix(c(X, K - X, M - X, N - M - K + X), nrow = 2, 
                   dimnames = list(c("anno", "notAnno"), c("group", 
                                                           "notGroup")))
    if (0) {
      p.value <- ifelse(all(cTab == 0), 1, stats::fisher.test(cTab, 
                                                              alternative = "greater")$p.value)
    }
    else {
      if (all(cTab == 0)) {
        p.value <- 1
      }
      else {
        if (p.tail == "one-tail") {
          p.value <- stats::fisher.test(cTab, alternative = "greater")$p.value
        }
        else {
          if (X >= K * M/N) {
            p.value <- stats::fisher.test(cTab, alternative = "greater")$p.value
          }
          else {
            p.value <- stats::fisher.test(cTab, alternative = "less")$p.value
          }
        }
      }
    }
    return(p.value)
  }
  doHypergeoTest <- function(genes.group, genes.term, genes.universe, 
                             p.tail) {
    genes.hit <- intersect(genes.group, genes.term)
    X <- length(genes.hit)
    K <- length(genes.group)
    M <- length(genes.term)
    N <- length(genes.universe)
    x <- X
    m <- M
    n <- N - M
    k <- K
    if (m == 0 || k == 0) {
      p.value <- 1
    }
    else {
      if (p.tail == "one-tail") {
        p.value <- stats::phyper(x, m, n, k, lower.tail = FALSE, 
                                 log.p = FALSE)
      }
      else {
        if (X >= K * M/N) {
          p.value <- stats::phyper(x, m, n, k, lower.tail = FALSE, 
                                   log.p = FALSE)
        }
        else {
          p.value <- stats::phyper(x, m, n, k, lower.tail = TRUE, 
                                   log.p = FALSE)
        }
      }
    }
    return(p.value)
  }
  doBinomialTest <- function(genes.group, genes.term, genes.universe, 
                             p.tail) {
    genes.hit <- intersect(genes.group, genes.term)
    X <- length(genes.hit)
    K <- length(genes.group)
    M <- length(genes.term)
    N <- length(genes.universe)
    if (K == 0 || M == 0 || N == 0) {
      p.value <- 1
    }
    else {
      if (p.tail == "one-tail") {
        p.value <- stats::pbinom(X, K, M/N, lower.tail = FALSE, 
                                 log.p = FALSE)
      }
      else {
        if (X >= K * M/N) {
          p.value <- stats::pbinom(X, K, M/N, lower.tail = FALSE, 
                                   log.p = FALSE)
        }
        else {
          p.value <- stats::pbinom(X, K, M/N, lower.tail = TRUE, 
                                   log.p = FALSE)
        }
      }
    }
    return(p.value)
  }
  zscoreHyper <- function(genes.group, genes.term, genes.universe) {
    genes.hit <- intersect(genes.group, genes.term)
    X <- length(genes.hit)
    K <- length(genes.group)
    M <- length(genes.term)
    N <- length(genes.universe)
    if (1) {
      x.exp <- K * M/N
      var.exp <- K * M/N * (N - M)/N * (N - K)/(N - 1)
      if (is.na(var.exp)) {
        z <- NA
      }
      else {
        if (var.exp != 0) {
          suppressWarnings(z <- (X - x.exp)/sqrt(var.exp))
        }
        else {
          z <- NA
        }
      }
    }
    else {
      x <- X
      m <- M
      n <- N - M
      k <- K
      suppressWarnings(d <- stats::dhyper(x, m, n, k, log = TRUE) - 
                         log(2))
      suppressWarnings(pupper <- stats::phyper(x, m, n, 
                                               k, lower.tail = FALSE, log.p = TRUE))
      suppressWarnings(plower <- stats::phyper(x - 1, m, 
                                               n, k, lower.tail = TRUE, log.p = TRUE))
      d[is.na(d)] <- -Inf
      pupper[is.na(pupper)] <- -Inf
      plower[is.na(plower)] <- -Inf
      a <- pupper
      b <- d - pupper
      a[b > 0] <- d[b > 0]
      b <- -abs(b)
      pmidupper <- a + log1p(exp(b))
      pmidupper[is.infinite(a)] <- a[is.infinite(a)]
      a <- plower
      b <- d - plower
      a[b > 0] <- d[b > 0]
      b <- -abs(b)
      pmidlower <- a + log1p(exp(b))
      pmidlower[is.infinite(a)] <- a[is.infinite(a)]
      up <- pmidupper < pmidlower
      if (any(up)) 
        z <- stats::qnorm(pmidupper, lower.tail = FALSE, 
                          log.p = TRUE)
      if (any(!up)) 
        z <- stats::qnorm(pmidlower, lower.tail = TRUE, 
                          log.p = TRUE)
    }
    return(z)
  }
  fcHyper <- function(genes.group, genes.term, genes.universe) {
    genes.hit <- intersect(genes.group, genes.term)
    X <- length(genes.hit)
    K <- length(genes.group)
    M <- length(genes.term)
    N <- length(genes.universe)
    x.exp <- K * M/N
    fc <- X/x.exp
    return(fc)
  }
  orFisher <- function(genes.group, genes.term, genes.universe) {
    genes.hit <- intersect(genes.group, genes.term)
    X <- length(genes.hit)
    K <- length(genes.group)
    M <- length(genes.term)
    N <- length(genes.universe)
    cTab <- matrix(c(X, K - X, M - X, N - M - K + X), nrow = 2, 
                   dimnames = list(c("anno", "notAnno"), c("group", 
                                                           "notGroup")))
    res <- stats::fisher.test(cTab)
    return(c(as.vector(res$estimate), as.vector(res$conf.int)))
  }
  if (verbose) {
    now <- Sys.time()
    message(sprintf("Next, prepare enrichment analysis (%s) ...", 
                    as.character(now)), appendLF = TRUE)
  }
  if (ontology.algorithm != "none") {
    background.annotatable.only <- TRUE
  }
  else {
    if (is.null(background.annotatable.only)) {
      if (length(background) == 0) {
        background.annotatable.only <- TRUE
      }
      else {
        background.annotatable.only <- FALSE
      }
    }
  }
  terms <- names(gs)
  if (background.annotatable.only) {
    genes.universe <- unique(unlist(gs[terms]))
  }
  else {
    genes.universe <- background
  }
  genes.group <- intersect(genes.universe, data)
  if (length(genes.group) == 0) {
    warnings("There is no gene being used.\n")
    return(NULL)
  }
  else {
    if (verbose) {
      now <- Sys.time()
      message(sprintf("\tThere are %d genes/SNPs of interest tested against %d genes/SNPs as the background (annotatable only? %s) (%s)", 
                      length(genes.group), length(genes.universe), 
                      background.annotatable.only, as.character(now)), 
              appendLF = TRUE)
    }
  }
  subg <- dnet::dDAGinduce(g = subg, nodes_query = terms, path.mode = path.mode)
  set_info <- data.frame(id = V(subg)$term_id, name = V(subg)$term_name, 
                         distance = V(subg)$term_distance, namespace = "term_namespace", 
                         row.names = V(subg)$name, stringsAsFactors = FALSE)
  if (ontology.algorithm == "none") {
    if (verbose) {
      now <- Sys.time()
      message(sprintf("Third, perform enrichment analysis using '%s' test (%s) ...", 
                      test, as.character(now)), appendLF = TRUE)
      if (is.null(which.distance)) {
        message(sprintf("\tThere are %d terms being used, each restricted within [%s] annotations", 
                        length(terms), paste(size.range, collapse = ",")), 
                appendLF = TRUE)
      }
      else {
        message(sprintf("\tThere are %d terms being used, each restricted within [%s] annotations and [%s] distance", 
                        length(terms), paste(size.range, collapse = ","), 
                        paste(which.distance, collapse = ",")), appendLF = TRUE)
      }
    }
    pvals <- sapply(terms, function(term) {
      genes.term <- unique(unlist(gs[term]))
      p.value <- switch(test, fisher = doFisherTest(genes.group, 
                                                    genes.term, genes.universe), hypergeo = doHypergeoTest(genes.group, 
                                                                                                           genes.term, genes.universe, p.tail = p.tail), 
                        binomial = doBinomialTest(genes.group, genes.term, 
                                                  genes.universe, p.tail = p.tail))
    })
    zscores <- sapply(terms, function(term) {
      genes.term <- unique(unlist(gs[term]))
      zscoreHyper(genes.group, genes.term, genes.universe)
    })
    fcs <- sapply(terms, function(term) {
      genes.term <- unique(unlist(gs[term]))
      fcHyper(genes.group, genes.term, genes.universe)
    })
    ls_or <- lapply(terms, function(term) {
      genes.term <- unique(unlist(gs[term]))
      orFisher(genes.group, genes.term, genes.universe)
    })
    df_or <- do.call(rbind, ls_or)
    ors <- df_or[, 1]
    CIl <- df_or[, 2]
    CIu <- df_or[, 3]
  }
  else if (ontology.algorithm == "pc" || ontology.algorithm == 
           "elim" || ontology.algorithm == "lea") {
    if (verbose) {
      now <- Sys.time()
      message(sprintf("Third, perform enrichment analysis based on '%s' test, and also using '%s' algorithm to respect ontology structure (%s) ...", 
                      test, ontology.algorithm, as.character(now)), 
              appendLF = TRUE)
    }
    if (verbose) {
      message(sprintf("\tThere are %d terms being used", 
                      length(V(subg))), appendLF = TRUE)
    }
    level2node <- dnet::dDAGlevel(subg, level.mode = "longest_path", 
                                  return.mode = "level2node")
    level2node.Hash <- list2env(level2node)
    nLevels <- length(level2node)
    node2pval.Hash <- new.env(hash = TRUE, parent = emptyenv())
    node2zscore.Hash <- new.env(hash = TRUE, parent = emptyenv())
    node2fc.Hash <- new.env(hash = TRUE, parent = emptyenv())
    node2or.Hash <- new.env(hash = TRUE, parent = emptyenv())
    node2CIl.Hash <- new.env(hash = TRUE, parent = emptyenv())
    node2CIu.Hash <- new.env(hash = TRUE, parent = emptyenv())
    if (ontology.algorithm == "pc") {
      for (i in nLevels:2) {
        currNodes <- get(as.character(i), envir = level2node.Hash, 
                         mode = "character")
        for (currNode in currNodes) {
          genes.term <- unique(unlist(gs[currNode]))
          pvalue_whole <- switch(test, fisher = doFisherTest(genes.group, 
                                                             genes.term, genes.universe), hypergeo = doHypergeoTest(genes.group, 
                                                                                                                    genes.term, genes.universe, p.tail = p.tail), 
                                 binomial = doBinomialTest(genes.group, genes.term, 
                                                           genes.universe, p.tail = p.tail))
          zscore_whole <- zscoreHyper(genes.group, genes.term, 
                                      genes.universe)
          fc_whole <- fcHyper(genes.group, genes.term, 
                              genes.universe)
          vec_whole <- orFisher(genes.group, genes.term, 
                                genes.universe)
          or_whole <- vec_whole[1]
          CIl_whole <- vec_whole[2]
          CIu_whole <- vec_whole[3]
          neighs.in <- igraph::neighborhood(subg, order = 1, 
                                            nodes = currNode, mode = "in")
          adjNodes <- setdiff(V(subg)[unlist(neighs.in)]$name, 
                              currNode)
          genes.parent <- unique(unlist(gs[adjNodes]))
          genes.group.parent <- intersect(genes.group, 
                                          genes.parent)
          genes.term.parent <- intersect(genes.term, 
                                         genes.parent)
          pvalue_relative <- switch(test, fisher = doFisherTest(genes.group.parent, 
                                                                genes.term.parent, genes.parent), hypergeo = doHypergeoTest(genes.group.parent, 
                                                                                                                            genes.term.parent, genes.parent, p.tail = p.tail), 
                                    binomial = doBinomialTest(genes.group.parent, 
                                                              genes.term.parent, genes.parent, p.tail = p.tail))
          zscore_relative <- zscoreHyper(genes.group.parent, 
                                         genes.term.parent, genes.parent)
          fc_relative <- fcHyper(genes.group.parent, 
                                 genes.term.parent, genes.parent)
          vec_relative <- orFisher(genes.group.parent, 
                                   genes.term.parent, genes.parent)
          or_relative <- vec_relative[1]
          CIl_relative <- vec_relative[2]
          CIu_relative <- vec_relative[3]
          pvalue <- max(pvalue_whole, pvalue_relative)
          assign(currNode, pvalue, envir = node2pval.Hash)
          zscore <- ifelse(pvalue_whole > pvalue_relative, 
                           zscore_whole, zscore_relative)
          assign(currNode, zscore, envir = node2zscore.Hash)
          fc <- ifelse(pvalue_whole > pvalue_relative, 
                       fc_whole, fc_relative)
          assign(currNode, fc, envir = node2fc.Hash)
          or <- ifelse(pvalue_whole > pvalue_relative, 
                       or_whole, or_relative)
          assign(currNode, or, envir = node2or.Hash)
          CIl <- ifelse(pvalue_whole > pvalue_relative, 
                        CIl_whole, CIl_relative)
          assign(currNode, CIl, envir = node2CIl.Hash)
          CIu <- ifelse(pvalue_whole > pvalue_relative, 
                        CIu_whole, CIu_relative)
          assign(currNode, CIu, envir = node2CIu.Hash)
        }
        if (verbose) {
          message(sprintf("\tAt level %d, there are %d nodes/terms", 
                          i, length(currNodes), appendLF = TRUE))
        }
      }
      root <- dnet::dDAGroot(subg)
      assign(root, 1, envir = node2pval.Hash)
      assign(root, 0, envir = node2zscore.Hash)
    }
    else if (ontology.algorithm == "elim") {
      sigNode2pval.Hash <- new.env(hash = TRUE, parent = emptyenv())
      ancNode2gene.Hash <- new.env(hash = TRUE, parent = emptyenv())
      if (is.null(elim.pvalue) || is.na(elim.pvalue) || 
          elim.pvalue > 1 || elim.pvalue < 0) {
        elim.pvalue <- 0.01
      }
      pval.cutoff <- elim.pvalue
      for (i in nLevels:1) {
        currNodes <- get(as.character(i), envir = level2node.Hash, 
                         mode = "character")
        currAnno <- gs[currNodes]
        for (currNode in currNodes) {
          genes.term <- unique(unlist(gs[currNode]))
          if (exists(currNode, envir = ancNode2gene.Hash, 
                     mode = "numeric")) {
            genes.elim <- get(currNode, envir = ancNode2gene.Hash, 
                              mode = "numeric")
            genes.term <- setdiff(genes.term, genes.elim)
          }
          pvalue <- switch(test, fisher = doFisherTest(genes.group, 
                                                       genes.term, genes.universe), hypergeo = doHypergeoTest(genes.group, 
                                                                                                              genes.term, genes.universe, p.tail = p.tail), 
                           binomial = doBinomialTest(genes.group, genes.term, 
                                                     genes.universe, p.tail = p.tail))
          zscore <- zscoreHyper(genes.group, genes.term, 
                                genes.universe)
          fc <- fcHyper(genes.group, genes.term, genes.universe)
          vec <- orFisher(genes.group, genes.term, genes.universe)
          or <- vec[1]
          CIl <- vec[2]
          CIu <- vec[3]
          assign(currNode, pvalue, envir = node2pval.Hash)
          assign(currNode, zscore, envir = node2zscore.Hash)
          assign(currNode, fc, envir = node2fc.Hash)
          assign(currNode, or, envir = node2or.Hash)
          assign(currNode, CIl, envir = node2CIl.Hash)
          assign(currNode, CIu, envir = node2CIu.Hash)
          if (pvalue < pval.cutoff) {
            assign(currNode, pvalue, envir = sigNode2pval.Hash)
            elimGenesID <- currAnno[[currNode]]
            dag.ancestors <- dnet::dDAGinduce(subg, currNode, 
                                              path.mode = "all_paths")
            ancestors <- setdiff(V(dag.ancestors)$name, 
                                 currNode)
            oldAncestors2GenesID <- sapply(ancestors, 
                                           function(ancestor) {
                                             if (exists(ancestor, envir = ancNode2gene.Hash, 
                                                        mode = "numeric")) {
                                               get(ancestor, envir = ancNode2gene.Hash, 
                                                   mode = "numeric")
                                             }
                                           })
            newAncestors2GenesID <- lapply(oldAncestors2GenesID, 
                                           function(oldGenes) {
                                             base::union(oldGenes, elimGenesID)
                                           })
            if (length(newAncestors2GenesID) > 0) {
              sapply(names(newAncestors2GenesID), function(ancestor) {
                assign(ancestor, newAncestors2GenesID[[ancestor]], 
                       envir = ancNode2gene.Hash)
              })
            }
          }
        }
        if (verbose) {
          num.signodes <- length(ls(sigNode2pval.Hash))
          num.ancnodes <- length(ls(ancNode2gene.Hash))
          num.elimgenes <- length(unique(unlist(as.list(ancNode2gene.Hash))))
          message(sprintf("\tAt level %d, there are %d nodes/terms: up to %d significant nodes, %d ancestral nodes changed (%d genes/SNPs eliminated)", 
                          i, length(currNodes), num.signodes, num.ancnodes, 
                          num.elimgenes), appendLF = TRUE)
        }
      }
    }
    else if (ontology.algorithm == "lea") {
      node2pvalo.Hash <- new.env(hash = TRUE, parent = emptyenv())
      if (is.null(lea.depth) || is.na(lea.depth) || lea.depth < 
          0) {
        lea.depth <- 2
      }
      depth.cutoff <- as.integer(lea.depth)
      for (i in nLevels:1) {
        currNodes <- get(as.character(i), envir = level2node.Hash, 
                         mode = "character")
        currAnno <- gs[currNodes]
        num.recalculate <- 0
        for (currNode in currNodes) {
          genes.term <- unique(unlist(gs[currNode]))
          pvalue.old <- switch(test, fisher = doFisherTest(genes.group, 
                                                           genes.term, genes.universe), hypergeo = doHypergeoTest(genes.group, 
                                                                                                                  genes.term, genes.universe, p.tail = p.tail), 
                               binomial = doBinomialTest(genes.group, genes.term, 
                                                         genes.universe, p.tail = p.tail))
          zscore.old <- zscoreHyper(genes.group, genes.term, 
                                    genes.universe)
          fc.old <- fcHyper(genes.group, genes.term, 
                            genes.universe)
          vec.old <- orFisher(genes.group, genes.term, 
                              genes.universe)
          or.old <- vec.old[1]
          CIl.old <- vec.old[2]
          CIu.old <- vec.old[3]
          assign(currNode, pvalue.old, envir = node2pvalo.Hash)
          neighs.out <- igraph::neighborhood(subg, order = depth.cutoff, 
                                             nodes = currNode, mode = "out")
          adjNodes <- setdiff(V(subg)[unlist(neighs.out)]$name, 
                              currNode)
          if (length(adjNodes) != 0) {
            if (1) {
              pvalue.children <- sapply(adjNodes, function(child) {
                if (exists(child, envir = node2pvalo.Hash, 
                           mode = "numeric")) {
                  get(child, envir = node2pvalo.Hash, 
                      mode = "numeric")
                }
              })
            }
            else {
              pvalue.children <- sapply(adjNodes, function(child) {
                if (exists(child, envir = node2pval.Hash, 
                           mode = "numeric")) {
                  get(child, envir = node2pval.Hash, 
                      mode = "numeric")
                }
              })
            }
            chNodes <- names(pvalue.children[pvalue.children < 
                                               pvalue.old])
            if (length(chNodes) > 0) {
              num.recalculate <- num.recalculate + 1
              genes.elim <- unique(unlist(gs[chNodes]))
              genes.term.new <- setdiff(genes.term, genes.elim)
              pvalue.new <- switch(test, fisher = doFisherTest(genes.group, 
                                                               genes.term.new, genes.universe), hypergeo = doHypergeoTest(genes.group, 
                                                                                                                          genes.term.new, genes.universe, p.tail = p.tail), 
                                   binomial = doBinomialTest(genes.group, 
                                                             genes.term.new, genes.universe, p.tail = p.tail))
              zscore.new <- zscoreHyper(genes.group, 
                                        genes.term.new, genes.universe)
              fc.new <- fcHyper(genes.group, genes.term.new, 
                                genes.universe)
              vec.new <- orFisher(genes.group, genes.term.new, 
                                  genes.universe)
              or.new <- vec.new[1]
              CIl.new <- vec.new[2]
              CIu.new <- vec.new[3]
              pvalue <- max(pvalue.new, pvalue.old)
              zscore <- ifelse(pvalue.new > pvalue.old, 
                               zscore.new, zscore.old)
              fc <- ifelse(pvalue.new > pvalue.old, fc.new, 
                           fc.old)
              or <- ifelse(pvalue.new > pvalue.old, or.new, 
                           or.old)
              CIl <- ifelse(pvalue.new > pvalue.old, 
                            CIl.new, CIl.old)
              CIu <- ifelse(pvalue.new > pvalue.old, 
                            CIu.new, CIu.old)
            }
            else {
              pvalue <- pvalue.old
              zscore <- zscore.old
              fc <- fc.old
              or <- or.old
              CIl <- CIl.old
              CIu <- CIu.old
            }
          }
          else {
            pvalue <- pvalue.old
            zscore <- zscore.old
            fc <- fc.old
            or <- or.old
            CIl <- CIl.old
            CIu <- CIu.old
          }
          assign(currNode, pvalue, envir = node2pval.Hash)
          assign(currNode, zscore, envir = node2zscore.Hash)
          assign(currNode, fc, envir = node2fc.Hash)
          assign(currNode, or, envir = node2or.Hash)
          assign(currNode, CIl, envir = node2CIl.Hash)
          assign(currNode, CIu, envir = node2CIu.Hash)
        }
        if (verbose) {
          message(sprintf("\tAt level %d, there are %d nodes/terms and %d being recalculated", 
                          i, length(currNodes), num.recalculate), appendLF = TRUE)
        }
      }
    }
    pvals <- unlist(as.list(node2pval.Hash))
    zscores <- unlist(as.list(node2zscore.Hash))
    fcs <- unlist(as.list(node2fc.Hash))
    ors <- unlist(as.list(node2or.Hash))
    CIl <- unlist(as.list(node2CIl.Hash))
    CIu <- unlist(as.list(node2CIu.Hash))
  }
  overlaps <- lapply(names(gs), function(term) {
    genes.term <- unique(unlist(gs[term]))
    x <- intersect(genes.group, genes.term)
    x
  })
  names(overlaps) <- names(gs)
  flag_filter <- sapply(overlaps, function(x) ifelse(length(x) >= 
                                                       min.overlap, TRUE, FALSE))
  if (sum(flag_filter) == 0) {
    warnings("It seems there are no terms meeting the specified 'size.range' and 'min.overlap'.\n")
    return(NULL)
  }
  gs <- gs[flag_filter]
  overlaps <- overlaps[flag_filter]
  common <- intersect(names(gs), names(zscores))
  ind_gs <- match(common, names(gs))
  ind_zscores <- match(common, names(zscores))
  gs <- gs[ind_gs[!is.na(ind_gs)]]
  overlaps <- overlaps[ind_gs[!is.na(ind_gs)]]
  zscores <- zscores[ind_zscores[!is.na(ind_zscores)]]
  fcs <- fcs[ind_zscores[!is.na(ind_zscores)]]
  pvals <- pvals[ind_zscores[!is.na(ind_zscores)]]
  ors <- ors[ind_zscores[!is.na(ind_zscores)]]
  CIl <- CIl[ind_zscores[!is.na(ind_zscores)]]
  CIu <- CIu[ind_zscores[!is.na(ind_zscores)]]
  flag <- !is.na(zscores)
  gs <- gs[flag]
  overlaps <- overlaps[flag]
  zscores <- zscores[flag]
  fcs <- fcs[flag]
  pvals <- pvals[flag]
  ors <- ors[flag]
  CIl <- CIl[flag]
  CIu <- CIu[flag]
  if (length(pvals) == 0) {
    warnings("There are no pvals being calculated.\n")
    return(NULL)
  }
  ind <- match(rownames(set_info), names(pvals))
  set_info <- set_info[!is.na(ind), ]
  zscores <- signif(zscores, digits = 3)
  fcs <- signif(fcs, digits = 3)
  pvals <- sapply(pvals, function(x) min(x, 1))
  ors <- signif(ors, digits = 3)
  CIl <- signif(CIl, digits = 3)
  CIu <- signif(CIu, digits = 3)
  if (verbose) {
    now <- Sys.time()
    message(sprintf("Last, adjust the p-values for %d terms (with %d minimum overlaps) using the %s method (%s) ...", 
                    length(pvals), min.overlap, p.adjust.method, as.character(now)), 
            appendLF = TRUE)
  }
  adjpvals <- stats::p.adjust(pvals, method = p.adjust.method)
  pvals <- signif(pvals, digits = 2)
  adjpvals <- sapply(adjpvals, function(x) min(x, 1))
  adjpvals <- signif(adjpvals, digits = 2)
  if (0) {
    tmp <- as.numeric(format(.Machine)["double.xmin"])
    tmp <- signif(tmp, digits = 2)
    pvals[pvals < tmp] <- tmp
    adjpvals[adjpvals < tmp] <- tmp
  }
  pvals <- sapply(pvals, function(x) {
    if (x < 0.1 & x != 0) {
      as.numeric(format(x, scientific = TRUE))
    }
    else {
      x
    }
  })
  adjpvals <- sapply(adjpvals, function(x) {
    if (x < 0.1 & x != 0) {
      as.numeric(format(x, scientific = TRUE))
    }
    else {
      x
    }
  })
  cross <- matrix(0, nrow = length(overlaps), ncol = length(overlaps))
  if (length(overlaps) >= 2) {
    for (i in seq(1, length(overlaps) - 1)) {
      x1 <- overlaps[[i]]
      for (j in seq(i + 1, length(overlaps))) {
        x2 <- overlaps[[j]]
        cross[i, j] <- length(intersect(x1, x2))
        cross[j, i] <- length(intersect(x1, x2))
      }
    }
    colnames(cross) <- rownames(cross) <- names(overlaps)
    diag(cross) <- sapply(overlaps, length)
  }
  eTerm <- list(term_info = set_info, annotation = gs, g = subg, 
                data = genes.group, background = genes.universe, overlap = overlaps, 
                fc = fcs, zscore = zscores, pvalue = pvals, adjp = adjpvals, 
                or = ors, CIl = CIl, CIu = CIu, cross = cross, call = match.call())
  class(eTerm) <- "eTerm"
  invisible(eTerm)
}







