rm(list = ls())
library(readr)
library(quanteda)
library(quanteda.textplots)
library(stm)
data <- read.csv("UNGDspeeches.csv")

# make a new col that has 1 if it is in depression years, 0 otherwise
depression_years <- c(1975, 1982, 1991, 2009)
data$depression <- ifelse(data$year%in%depression_years, 1, 0)

# preprocessing data (general, nd(not depression), d(depression year))
quant_corpus <- corpus(data, text_field = "text")
quant_corpus_nd <- corpus(data[data$depression == 0,], text_field = "text")
quant_corpus_d <- corpus(data[data$depression == 1,], text_field = "text")

# tokenize the corpus for all years, non depression years, depression years
toks <- tokens(quant_corpus, remove_punct = TRUE, remove_numbers = TRUE)
toks_nd <- tokens(quant_corpus_nd, remove_punct = TRUE, remove_numbers = TRUE)
toks_d <- tokens(quant_corpus_d, remove_punct = TRUE, remove_numbers = TRUE)

# stemming words
toks <- tokens_wordstem(toks)
toks_nd <- tokens_wordstem(toks_nd)
toks_d <- tokens_wordstem(toks_d)

# remove stop words
toks <- tokens_select(toks,  stopwords("en"), selection = "remove")
toks_nd <- tokens_select(toks_nd,  stopwords("en"), selection = "remove")
toks_d <- tokens_select(toks_d,  stopwords("en"), selection = "remove")

# create dfm
dfm <- dfm(toks)
dfm_nd <- dfm(toks_nd)
dfm_d <- dfm(toks_d)

# trim dfm to only include words in min .05 of docs
dfm_trimmed <- dfm_trim(dfm, min_docfreq = .05, docfreq_type = "prop")
dfm_trimmed_nd <- dfm_trim(dfm_nd, min_docfreq = .05, docfreq_type = "prop")
dfm_trimmed_d <- dfm_trim(dfm_d, min_docfreq = .05, docfreq_type = "prop")


# word cloud
textplot_wordcloud(dfm_trimmed, col = "black")
textplot_wordcloud(dfm_trimmed_nd, col = "black")
textplot_wordcloud(dfm_trimmed_d, col = "black")


# fightin words
clusterFightinWords <- function(dfm, clust.vect, alpha.0=100) {
  # we need to get the overall corpus word distribution and the cluster-specific words dists
  # y_{kw} in Monroe et al. 
  overall.terms <- colSums(dfm)
  # n and n_k in Monroe et al. 
  n <- sum(overall.terms)
  # alpha_{kw} in Monroe et al. 
  prior.terms <- overall.terms / n * alpha.0
  # y_{kw}(i) in Monroe et al.
  cluster.terms <- colSums(dfm[clust.vect, ])
  # n_k(i) in Monroe et al.
  cluster.n <- sum(cluster.terms)
  
  cluster.term.odds <- 
    (cluster.terms + prior.terms) / 
    (cluster.n + alpha.0 - cluster.terms - prior.terms)
  overall.term.odds <- 
    (overall.terms + prior.terms) / 
    (n + alpha.0 - overall.terms - prior.terms)
  
  log.odds <- log(cluster.term.odds) - log(overall.term.odds)
  
  variance <- 1/(cluster.terms + prior.terms) + 1/(overall.terms + prior.terms)
  
  # return the variance weighted log-odds for each term
  output <- log.odds / sqrt(variance)
  names(output) <- colnames(dfm)
  return(output)
}

# find most defining words (non depression words)
nd_words <- clusterFightinWords(dfm_trimmed, data$depression == 0)
sort(nd_words, decreasing = T)[5:25]

# (depression words)
d_words <- clusterFightinWords(dfm_trimmed, data$depression == 1)
sort(d_words, decreasing = T)[1:20]




# preprocess
temp <- textProcessor(documents = data$text, metadata = data)
out <- prepDocuments(temp$documents, temp$vocab, temp$meta)

# build model 
model.stm <- stm(out$documents, out$vocab, K = 10, prevalence = ~depression, 
                 data = out$meta, max.em.its = 10)


# get topics
labelTopics(model.stm)

# find important thoughts
findThoughts(model.stm, texts = out$meta$text, topics = 10, n = 3)

# find which topics belong to depresssion/nondepression years
model.stm.ee <- estimateEffect(1:10~depression, model.stm, meta = out$meta)

# plot
plot(model.stm.ee, "depression", method = "difference", cov.value1 = 1, cov.value2 = 0)
