#
# HMM model
#
library(rqhmm)

# NOTE: I would like to do the HMM with a 'scaled' neg. binomial
#       but I have not implemented that yet, so first run will 
#       use a DGamma instead
#


# 1. no RNA info (use TSS prior)
#
# state paths:
#    B -> T -> PP -> D -> B
#           -> B [optional link]
#
#  state  #pro  dist
#    B     [1]   [1]
#    T     [2]   [2]
#   PP     [3]   [2]
#    D     [4]   [2]
#
# #pro: [1] background level (Poisson)
#       [2] body level (NegBinom | DGamma)
#       [3] fixed scale up of body level (NegBinom | DGamma)
#       [4] fixed scale down of body level (NegBinom | DGamma)
#
# dist: [1] background distances (Geom0)
#       [2] body distances (Geom0)
#
splithmm1.hmm <- function(scaleTP, scalePD, with.shortcut = FALSE, use.negbinom = FALSE) {
  N = 4
  
  # B, T, P, D
  
  # valid transitions
  vtbl = NULL
  tnames = c("autocorr_covar", "autocorr", "autocorr", "autocorr")
  
  if (!with.shortcut) {
    #
    # no link from T to B directly
    #
    vtbl = rbind(
      c(1, 2, 0, 0), # B to: B or T
      c(0, 1, 2, 0), # T to: T or P
      c(0, 0, 1, 2), # P to: P or D
      c(2, 0, 0, 1)) # D to: D or B
  } else {
    #
    # no link from T to B directly
    #
    vtbl = rbind(
      c(1, 2, 0, 0), # B to: B or T
      c(3, 1, 2, 0), # T to: T or B or P
      c(0, 0, 1, 2), # P to: P or D
      c(2, 0, 0, 1)) # D to: D or B
  }

  tdist = "dgamma"
  if (use.negbinom)
    tdist = "neg_binomial"

  # transition groups
  tgrps = list(2:N)  # share size distribution; it's a poor fit for all, but want to avoid getting a tiny 'T' state
  if (with.shortcut)
   tgrps = list(3:N) # must exclude T=2 state as it has a different number of transitions

  # emission groups
  egrps = new.emission.groups(N, 2)
  egrps = add.emission.groups(egrps, states = c(2, 3, 4), slots = c(2, 2, 2)) # share GROseq over T, P and D [scaled]

  hmm = new.qhmm(list(c(1, 1), 1), # covar: TSS signal
    vtbl,
    tnames,
    list(c("geometric", "poisson"),  # B
         c("geometric", tdist),      # T
         c("geometric", tdist),      # P
         c("geometric", tdist)),     # D
    transition.groups = tgrps,
    emission.groups = egrps)
  
  # set initial parameters
  set.initial.probs.qhmm(hmm, c(1, 0, 0, 0)) # start with background
  
  # set transitions
  set.transition.params.qhmm(hmm, 2:N, 0.99)
  
  # set emissions
  geom.base = 0
  # background
  set.emission.params.qhmm(hmm, 1, 1/3000, slot = 1)
  set.emission.option.qhmm(hmm, 1, "base", geom.base, slot = 1)
  set.emission.params.qhmm(hmm, 1, 0.1, slot = 2) # lambda
  
  # transcribed
  for (i in 2:N) {
    set.emission.params.qhmm(hmm, i, 1/10, slot = 1)
    set.emission.option.qhmm(hmm, i, "base", geom.base, slot = 1)
  }
  
  if (use.negbinom) {
    set.emission.params.qhmm(hmm, 2:N, c(10, 0.1), slot = 2)
    # set scale factors
    stop("not implemented!")
  } else {
    set.emission.params.qhmm(hmm, 2:N, c(1, 1), slot = 2) # gamma params
    set.emission.option.qhmm(hmm, 3, "scale_private", scaleTP, slot = 2)
    set.emission.option.qhmm(hmm, 4, "scale_private", scalePD, slot = 2)
  }

  return(hmm)
}


# 2. w/ RNA info (use TSS prior)
#    B -> U
#      -> FEXON -> INTRON -> EXON -> INTRON
#                                 -> PPAUSE -> PDECAY -> B
#
#   state  #pro  dist  #rna
#     B     [1]   [1]   [1]
#     U     [2]   [2]   [1]
#   FEXON   [2]   [2]   [2]
#  INTRON   [3]   [2]   [1]
#   EXON    [3]   [2]   [2]
#  PPAUSE   [4]   [2]   [1]
#  PDECAY   [5]   [2]   [1]
#
# #pro: [1] background level (Poisson)
#       [2] start level (NegBinom | DGamma)
#       [3] body level (NegBinom | DGamma)
#       [4] fixed scale up of body level (NegBinom | DGamma)
#       [5] fixed scale down of body level (NegBinom | DGamma)
#
# dist: [1] background distances (Geom0)
#       [2] body distances (Geom0)
#
# #rna: [1] background level (Poisson)
#       [2] exon level (NegBinom | DGamma)
#
splithmm2.hmm <- function(scaleTP, scalePD, use.negbinom = FALSE) {
  N = 7
  
  # B, U, F, I, E, P, D
  
  # valid transitions
  vtbl = NULL
  tnames = c("autocorr_covar", rep("autocorr", N - 1))

  vtbl = rbind(
    c(1, 2, 3, 0, 0, 0, 0), # B to: B or U or F
    c(2, 1, 0, 0, 0, 0, 0), # U to: U or B
    c(0, 0, 1, 2, 0, 0, 0), # F to: F or I
    c(0, 0, 0, 1, 2, 0, 0), # I to: I or E
    c(0, 0, 0, 2, 1, 3, 0), # E to: E or I or P
    c(0, 0, 0, 0, 0, 1, 2), # P to: P or D
    c(2, 0, 0, 0, 0, 0, 1)) # D to: D or B

  tdist = "dgamma"
  if (use.negbinom)
    tdist = "neg_binomial"

  # emission groups
  egrps = new.emission.groups(N, 3)
  # distances
  egrps = add.emission.groups(egrps, states = 2:7, slots = rep(1, 6)) # share distance distribution between all non-background states
  # pro/gro
  egrps = add.emission.groups(egrps, states = 4:7, slots = rep(2, 4)) # share GROseq over E, I, P, D (last two scaled versions)
  egrps = add.emission.groups(egrps, states = 2:3, slots = c(2, 2)) # share GROseq over U and F
  # rna
  egrps = add.emission.groups(egrps, states = c(1, 2, 4, 6, 7), slots = c(3, 3, 3, 3, 3)) # share background RNAseq (B, U, I, P, D)
  egrps = add.emission.groups(egrps, states = c(3, 5), slots = c(3, 3)) # share exon RNAseq (F, E)

  hmm = new.qhmm(list(c(1, 1, 1), 1), # covar: TSS signal
    vtbl,
    tnames,
    list(c("geometric", "poisson", "poisson"),  # B
         c("geometric", tdist, "poisson"),      # U
         c("geometric", tdist, tdist),          # F
         c("geometric", tdist, "poisson"),      # I
         c("geometric", tdist, tdist),          # E
         c("geometric", tdist, "poisson"),      # P
         c("geometric", tdist, "poisson")),     # D
    transition.groups = list(c(3, 5)), # share size distribution between F and E (both exons)
    emission.groups = egrps)
  
  # set initial parameters
  set.initial.probs.qhmm(hmm, c(1, rep(0, N - 1))) # start with background
  
  # set transitions
  set.transition.params.qhmm(hmm, 2:N, 0.99)

  # set emissions
  geom.base = 0
  # background
  set.emission.params.qhmm(hmm, 1, 1/3000, slot = 1)
  set.emission.option.qhmm(hmm, 1, "base", geom.base, slot = 1)
  set.emission.params.qhmm(hmm, 1, 0.1, slot = 2) # lambda
  set.emission.params.qhmm(hmm, 1, 0.1, slot = 3) # lambda
  
  # transcribed
  
  # dist
  for (i in 2:N) {
    set.emission.params.qhmm(hmm, i, 1/10, slot = 1)
    set.emission.option.qhmm(hmm, i, "base", geom.base, slot = 1)
  }
  
  # pro/gro
  if (use.negbinom) {
    set.emission.params.qhmm(hmm, 2:N, c(10, 0.1), slot = 2)
    # set scale factors
    stop("not implemented!")
  } else {
    set.emission.params.qhmm(hmm, 2:N, c(1, 1), slot = 2) # gamma params
    set.emission.option.qhmm(hmm, 6, "scale_private", scaleTP, slot = 2)
    set.emission.option.qhmm(hmm, 7, "scale_private", scalePD, slot = 2)
  }
  
  # rna
  set.emission.params.qhmm(hmm, c(2,4,6,7) , 0.1, slot = 3) # lambda
  
  if (use.negbinom) {
    set.emission.params.qhmm(hmm, c(3, 5), c(10, 0.1), slot = 3)
  } else {
    set.emission.params.qhmm(hmm, c(3, 5), c(1, 1), slot = 3) # gamma params
  }

  return(hmm)
}

#
# 3. three state HMM
#
# B -> T -> D
#
# state paths:
#    B -> T -> D -> B
#           -> B (optional link)
#
#  state  #pro  dist
#    B     [1]   [1]
#    T     [2]   [2]
#    D     [3]   [2]
#
# #pro: [1] background level (Poisson)
#       [2] body level (NegBinom | DGamma)
#       [3] fixed scale down of body level (NegBinom | DGamma)
#
# dist: [1] background distances (Geom0)
#       [2] body distances (Geom0)
#
splithmm3.hmm <- function(scaleTD, with.shortcut = FALSE, no.egrps = FALSE, poisson.decay = FALSE, use.negbinom = FALSE) {
  N = 3
  
  if (poisson.decay & !no.egrps)
    stop("poisson.decay requires no.egrps")
  
  # B, T, D
  
  # valid transitions
  vtbl = NULL
  tnames = c("autocorr_covar", "autocorr", "autocorr")
  
  if (!with.shortcut) {
    #
    # no link from T to B directly
    #
    vtbl = rbind(
      c(1, 2, 0), # B to: B or T
      c(0, 1, 2), # T to: T or D
      c(2, 0, 1)) # D to: D or B
  } else {
    #
    # no link from T to B directly
    #
    vtbl = rbind(
      c(1, 2, 0), # B to: B or T
      c(3, 1, 2), # T to: T or D or B
      c(2, 0, 1)) # D to: D or B
  }

  tdist = "dgamma"
  if (use.negbinom)
    tdist = "neg_binomial"
  ddist = tdist
  if (poisson.decay)
    ddist = "poisson"

  # transition groups
  tgrps = list(2:N)  # share size distribution; it's a poor fit for all, but want to avoid getting a tiny 'T' state
  if (with.shortcut)
   tgrps = NULL # must exclude T=2 state as it has a different number of transitions

  # emission groups
  egrps = NULL
  if (!no.egrps) {
    egrps = new.emission.groups(N, 2)
    egrps = add.emission.groups(egrps, states = c(2, 3), slots = c(2, 2)) # share GROseq over T and D [scaled]
    egrps = add.emission.groups(egrps, states = c(2, 3), slots = c(1, 1)) # share distance over T and D
  }

  hmm = new.qhmm(list(c(1, 1), 1), # covar: TSS signal
    vtbl,
    tnames,
    list(c("geometric", "poisson"),  # B
         c("geometric", tdist),      # T
         c("geometric", ddist)),     # D
    transition.groups = tgrps,
    emission.groups = egrps)
  
  # set initial parameters
  set.initial.probs.qhmm(hmm, c(1, 0, 0)) # start with background
  
  # set transitions
  set.transition.params.qhmm(hmm, 2:N, 0.99)
  
  # set emissions
  geom.base = 0
  # background
  set.emission.params.qhmm(hmm, 1, 1/3000, slot = 1)
  set.emission.option.qhmm(hmm, 1, "base", geom.base, slot = 1)
  set.emission.params.qhmm(hmm, 1, 0.1, slot = 2) # lambda
  
  # transcribed
  for (i in 2:N) {
    set.emission.params.qhmm(hmm, i, 1/10, slot = 1)
    set.emission.option.qhmm(hmm, i, "base", geom.base, slot = 1)
  }
  
  if (use.negbinom) {
    if (!no.egrps) {
      set.emission.params.qhmm(hmm, 2:N, c(10, 0.1), slot = 2)
      # set scale factors
      stop("not implemented!")
    } else {
      if (poisson.decay) {
        set.emission.params.qhmm(hmm, 2, c(10, 0.1), slot = 2)
        set.emission.params.qhmm(hmm, 3, 0.1, slot = 2)
      } else {
        set.emission.params.qhmm(hmm, 2, c(10, 0.1), slot = 2)
        set.emission.params.qhmm(hmm, 2, c(0.1, 0.1), slot = 2)
      }
    }
  } else {
    if (!no.egrps) {
      set.emission.params.qhmm(hmm, 2:N, c(1, 1), slot = 2) # gamma params
      set.emission.option.qhmm(hmm, 3, "scale_private", scaleTD, slot = 2)
    } else {
      if (poisson.decay) {
        set.emission.params.qhmm(hmm, 2, c(1, 1), slot = 2) # gamma params
        set.emission.params.qhmm(hmm, 3, 0.1, slot = 2) # gamma params
      } else {
        set.emission.params.qhmm(hmm, 2, c(1, 1), slot = 2) # gamma params
        set.emission.params.qhmm(hmm, 3, c(1, 0.1), slot = 2) # gamma params 
      }
    }
  }

  return(hmm)
}

#
# 4. five state HMM
#
# B -> I -> P -> T -> D
#
# state paths:
#    B -> I -> P -> T -> D -> B
#           -....-> T -> D -> B [skippable pause]
#                     -> B (optional link)
#
#  state  #pro  dist
#    B     [1]   [1]
#    I     [2]   [2]
#    P     [3]   [2]
#    T     [4]   [3]
#    D     [5]   [4]
#
# #pro: [1] background level (Poisson)
#       [2] initiation level (Poisson)
#       [4] pause level (Poisson)
#       [4] body level (NegBinom | DGamma)
#       [5] decay level (Poisson)
#
# dist: [1] background distances (Geom0)
#       [2] initiation/pause distances (Geom0)
#       [3] body distances (Geom0)
#       [4] decay distances (Geom0)
#
splithmm4.hmm <- function(with.shortcut = TRUE, no.egrps = FALSE, use.negbinom = FALSE) {
  N = 5

  # B, I, P, T, D
  
  # valid transitions
  vtbl = NULL
  tnames = c("autocorr_covar", rep("autocorr", N - 1))
  
  if (!with.shortcut) {
    #
    # no link from T to B directly
    #
    vtbl = rbind(
      c(1, 2, 0, 0, 0), # B to: B or I
      c(0, 1, 2, 3, 0), # I to: I or P or T
      c(0, 0, 1, 2, 0), # P to: P or T
      c(0, 0, 0, 1, 2), # T to: T or D
      c(2, 0, 0, 0, 1)) # D to: D or B
  } else {
    #
    # no link from T to B directly
    #
    vtbl = rbind(
      c(1, 2, 0, 0, 0), # B to: B or I
      c(0, 1, 2, 3, 0), # I to: I or P or T
      c(0, 0, 1, 2, 0), # P to: P or T
      c(3, 0, 0, 1, 2), # T to: T or D or B
      c(2, 0, 0, 0, 1)) # D to: D or B
  }
  
  tdist = "dgamma"
  if (use.negbinom)
    tdist = "neg_binomial"
  ddist = "poisson"

  # transition groups
  
# NOT valid!!
#  tgrps = list(2:3, 4:5)  # share size distribution; it's a poor fit for all, but want to avoid getting a tiny 'T' state
#  if (with.shortcut)
#   tgrps = list(2:3) # must exclude T=4 state as it has a different number of transitions
  tgrps = list(4:5)
  if (with.shortcut)
    tgrps = NULL # must exclude T=4 state as it has a different number of transitions

  # emission groups
  egrps = NULL
  if (!no.egrps) {
    egrps = new.emission.groups(N, 2)
    egrps = add.emission.groups(egrps, states = c(2, 3), slots = c(1, 1)) # share distance over I and P
    egrps = add.emission.groups(egrps, states = c(2, 5), slots = c(2, 2)) # share reads over I and D
  }

  hmm = new.qhmm(list(c(1, 1), 1), # covar: TSS signal
    vtbl,
    tnames,
    list(c("geometric", "poisson"),  # B
         c("geometric", "poisson"),  # I
         c("geometric", "poisson"),  # P
         c("geometric", tdist),      # T
         c("geometric", ddist)),     # D
    transition.groups = tgrps,
    emission.groups = egrps)
  
  # set initial parameters
  set.initial.probs.qhmm(hmm, c(1, rep(0, N - 1))) # start with background

  # set transitions
  set.transition.params.qhmm(hmm, 2:3, 0.50)
  set.transition.params.qhmm(hmm, 4:N, 0.99)
  
  # set emissions
  geom.base = 0
  # background
  set.emission.params.qhmm(hmm, 1, 1/3000, slot = 1)
  set.emission.option.qhmm(hmm, 1, "base", geom.base, slot = 1)
  set.emission.params.qhmm(hmm, 1, 0.1, slot = 2) # lambda
  
  # initiation
  set.emission.params.qhmm(hmm, 2, 1/10, slot = 1)
  set.emission.option.qhmm(hmm, 2, "base", geom.base, slot = 1)
  set.emission.params.qhmm(hmm, 2, 0.2, slot = 2) # lambda
  
  # pause
  set.emission.params.qhmm(hmm, 3, 1/10, slot = 1)
  set.emission.option.qhmm(hmm, 3, "base", geom.base, slot = 1)
  set.emission.params.qhmm(hmm, 3, 10, slot = 2) # lambda
  
  # transcribed
  for (i in 4:N) {
    set.emission.params.qhmm(hmm, i, 1/10, slot = 1)
    set.emission.option.qhmm(hmm, i, "base", geom.base, slot = 1)
  }
  
  if (use.negbinom) {
    set.emission.params.qhmm(hmm, 4, c(10, 0.1), slot = 2)
    set.emission.params.qhmm(hmm, 5, 0.1, slot = 2)
  } else {
    set.emission.params.qhmm(hmm, 4, c(1, 1), slot = 2) # gamma params
    set.emission.params.qhmm(hmm, 5, 0.1, slot = 2) # gamma params
  }

  return(hmm)
}
