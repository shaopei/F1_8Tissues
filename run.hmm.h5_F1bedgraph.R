#
# New HMM from work on dense k562 work
#
#

scriptPath <- function() {
        cmdArgs <- commandArgs(trailingOnly = FALSE)
        needle <- "--file="
        match <- grep(needle, cmdArgs)
        if (length(match) > 0) {
                # Rscript
                return(normalizePath(dirname(sub(needle, "", cmdArgs[match]))))
        } else {
                # 'source'd via R console
                return(normalizePath(dirname(sys.frames()[[1]]$ofile)))
        }
}

library(tunits)
source(paste(scriptPath(), "hmm.prototypes.R", sep="/"))

step = 50

args <- commandArgs(trailingOnly=TRUE)

name = args[1]

bed.path = args[2]
bwPlus.path = args[3]
bwMinus.path = args[4]
ref.params.path = args[5] #optional


chrom=c("chr1",  "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17","chr18", "chr19", 
  "chr2",  "chr3",  "chr4" , "chr5" , "chr6" , "chr7" , "chr8", "chr9", "chrX")
#
# prepare data
#
dregBED.covar <- function(dataset, dreg.bed) {
  N = length(dataset)
  covars = vector(mode="list", length=N)
  
  for (i in 1:N) {
    chrom = dataset[[i]]$chrom
    positions = dataset[[i]]$positions
    step = dataset[[i]]$step

    bed.i = dreg.bed[dreg.bed[,1] == chrom,]

    covars.i = 1 - sapply(positions, function(pos) {
      pos = (pos - 1) * step
      idx = which(bed.i[,2] <= pos & bed.i[,3] > pos)

      if (length(idx) >= 1)
        return(max(bed.i[idx, 4])) # use 4 for bedgraph. 5 for 6 bed
      return(0)
    })

    covars[[i]] = covars.i
  }

  return(covars)
}

dreg.clamp <- function(lst) lapply(lst, function(values) {
  1 - pmax(0, pmin(1 - values, 1))
})



all.dset = load.dataset(chrom, list(bw.plus = bwPlus.path, bw.minus = bwMinus.path), step, transform=TRUE)
train.dset = train.dataset(all.dset)

covars = dregBED.covar(all.dset, read.table(bed.path))
covars.clamp = dreg.clamp(covars)

# extended
all.dset.c <- all.dset
for (i in 1:length(all.dset)) {
  all.dset.c[[i]]$covar = covars.clamp[[i]]
}

#
# create HMM instance
hmm.dreg = splithmm4.hmm(with.shortcut = TRUE, no.egrps = FALSE, use.negbinom = TRUE)
#
# if present, pre-load starting parameters
if (!is.na(ref.params.path)) {
  load(ref.params.path) # defines hmm.params
  restore.params.qhmm(hmm.dreg, hmm.params)
}

# run
trace.dreg = em.qhmm(hmm.dreg, train.dset, covar.lst = covars.clamp, n_threads = 2)
if (length(trace.dreg$loglik) == 2) {
  # give it another change in case some error occured (I should fix this on the QHMM side ...
  trace.dreg = em.qhmm(hmm.dreg, train.dset, covar.lst = covars.clamp, n_threads = 2)
}


# save parameters
hmm.params = collect.params.qhmm(hmm.dreg)
save(hmm.params, file=paste(name, "_h5.params.Rdata", sep=''))

# save predictions

# 1. just body
preds.dreg = decode.dataset(hmm.dreg, all.dset.c, 2, covar.lst = covars.clamp)
write.track(preds.dreg, name, name, paste(name, "_h5.preds.bed", sep=''))

# 2. full version
preds.full.dreg = decode.dataset(hmm.dreg, all.dset.c, 2:5, covar.lst = covars.clamp)
write.track(preds.full.dreg, name, name, paste(name, "_h5.preds.full.bed", sep=''))

# 3. extended version
write.extended.track(preds.full.dreg, preds.dreg, name, name, paste(name, "_h5.preds.ext.bed", sep=''))
