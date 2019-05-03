### ABBA-BABA test for Cyrtandra data ###

require(HybridCheck)

# Kauai
Analysis1 <- HC$new("Kauai-concat.NOGAPS.fa")

populations1 <- list(
  P1 = c('Clong'),
  P2 = c('ClXk'),
  P3 = c('Ckaua'),
  Others = c('CaffW','Cwawr','Cwawr2')
)

Analysis1$setPopulations(populations1)

Analysis1$prepareFourTaxonTests()

popCombos <- list(
  c(P1 = "P3", P2 = "P2", P3 = "P1", A = "Others")
)

Analysis1$prepareFourTaxonTests(popCombos)

Analysis1$runFourTaxonTests(selections = 'ALL', 
                            blockLength = 1000L)

fttResults1 <- Analysis1$tabulateFourTaxonTests(selections='TESTED',
                                                neat=T,
                                                global=T)

split_list <- split(fttResults1, fttResults1$P1)
set1 <- as.data.frame(split_list[[1]])
set2 <- as.data.frame(split_list[[2]])

# set1 = Ckaua, ClXk, Clong, Others
# set2 = Clong, ClXk, Ckaua, Others

set1$ABBA   # 425
set1$BABA   # 281

set2$ABBA  # 840
set2$BABA  # 281

set1$X2_P # 2.6e-04
set2$X2_P # 5.4e-38

set1$D_jEstimate  # 0.204
set2$D_jEstimate  # 0.499

set1$Fd_1DD4_jEstimate  # 0.026
set1$Fd_D2D4_jEstimate  # 0

set2$Fd_1DD4_jEstimate  # 0.099
set2$Fd_D2D4_jEstimate  # 0

################################
# Oahu 1
Analysis1 <- HC$new("Oahu-concat.NOGAPS.fa")

populations1 <- list(
  P1 = c('CcXh'),
  P2 = c('Cgran','Ckaul','Chaw3','Ccalp'),
  P3 = c('CsXg'),
  Others = c('Ccord')
)

Analysis1$setPopulations(populations1)

Analysis1$prepareFourTaxonTests()

popCombos <- list(
  c(P1 = "P3", P2 = "P2", P3 = "P1", A = "Others")
)

Analysis1$prepareFourTaxonTests(popCombos)

Analysis1$runFourTaxonTests(selections = 'ALL', 
                            blockLength = 1000L)

fttResults1 <- Analysis1$tabulateFourTaxonTests(selections='TESTED',
                                                neat=T,
                                                global=T)

split_list <- split(fttResults1, fttResults1$P2)
set1 <- as.data.frame(split_list[[1]])
set2 <- as.data.frame(split_list[[2]])

# set1 = P3, Others, P1, P2
# set2 = P1, P2, P3, Others

set1$ABBA   # 376
set1$BABA   # 300

set2$ABBA  # 697
set2$BABA  # 651

set1$X2_P # 8.3e-09
set2$X2_P # 3.7e-114

set1$D_jEstimate  # 0.113
set2$D_jEstimate  # 0.034

set1$Fd_1DD4_jEstimate  # 0.014
set1$Fd_D2D4_jEstimate  # 0

set2$Fd_1DD4_jEstimate  # 0.009
set2$Fd_D2D4_jEstimate  # 0

#############################
### Oahu 2
Analysis1 <- HC$new("Oahu-concat.NOGAPS.fa")

populations1 <- list(
  P1 = c('Cgran'),
  P2 = c('Ckaul'),
  P3 = c('Chaw3'),
  Others = c('Ccalp')
)

Analysis1$setPopulations(populations1)

Analysis1$prepareFourTaxonTests()

popCombos <- list(
  c(P1 = "P3", P2 = "P2", P3 = "P1", A = "Others")
)

Analysis1$prepareFourTaxonTests(popCombos)

Analysis1$runFourTaxonTests(selections = 'ALL', 
                            blockLength = 1000L)

fttResults1 <- Analysis1$tabulateFourTaxonTests(selections='TESTED',
                                                neat=T,
                                                global=T)

split_list <- split(fttResults1, fttResults1$P3)
set1 <- as.data.frame(split_list[[1]])
set2 <- as.data.frame(split_list[[2]])

# set1 = P1, P2, Others, P3
# set2 = P1, P2, P3, Others

set1$ABBA   # 774
set1$BABA   # 814

set2$ABBA  # 814
set2$BABA  # 774

set1$X2_P # 5.2e-46
set2$X2_P # 5.2e-46

set1$D_jEstimate  # -0.025
set2$D_jEstimate  # 0.025

set1$Fd_1DD4_jEstimate  # 0
set1$Fd_D2D4_jEstimate  # 0.008

set2$Fd_1DD4_jEstimate  # 0.001
set2$Fd_D2D4_jEstimate  # 0

#############################
### Maui Nui 1
Analysis1 <- HC$new("Maui-Nui-concat.NOGAPS.fa")

populations1 <- list(
  P1 = c('Chaw2'),
  P2 = c('Chaw1'),
  P3 = c('Cpalu'),
  Others = c('Cplat','Cmunr','ChXg')
)

Analysis1$setPopulations(populations1)

Analysis1$prepareFourTaxonTests()

popCombos <- list(
  c(P1 = "P3", P2 = "P2", P3 = "P1", A = "Others")
)

Analysis1$prepareFourTaxonTests(popCombos)

Analysis1$runFourTaxonTests(selections = 'ALL', 
                            blockLength = 1000L)

fttResults1 <- Analysis1$tabulateFourTaxonTests(selections='TESTED',
                                                neat=T,
                                                global=T)

split_list <- split(fttResults1, fttResults1$P3)
set1 <- as.data.frame(split_list[[1]])
set2 <- as.data.frame(split_list[[2]])

# set1 = P1, P2, P3, Others
# set2 = P3, P2, P1, Others

set1$ABBA   # 413
set1$BABA   # 425

set2$ABBA  # 249
set2$BABA  # 202

set1$X2_P # 1.5e-80
set2$X2_P # 3.5e-09

set1$D_jEstimate  # -0.015
set2$D_jEstimate  # 0.103

set1$Fd_1DD4_jEstimate  # 0
set1$Fd_D2D4_jEstimate  # 0.018

set2$Fd_1DD4_jEstimate  # 0.012
set2$Fd_D2D4_jEstimate  # 0

### Maui Nui 2
Analysis1 <- HC$new("Maui-Nui-concat.NOGAPS.fa")

populations1 <- list(
  P1 = c('Chaw2'),
  P2 = c('Chaw1'),
  P3 = c('Cpalu'),
  Others = c('Cplat','Cmunr')
)

Analysis1$setPopulations(populations1)

Analysis1$prepareFourTaxonTests()

popCombos <- list(
  c(P1 = "P3", P2 = "P2", P3 = "P1", A = "Others")
)

Analysis1$prepareFourTaxonTests(popCombos)

Analysis1$runFourTaxonTests(selections = 'ALL', 
                            blockLength = 1000L)

fttResults1 <- Analysis1$tabulateFourTaxonTests(selections='TESTED',
                                                neat=T,
                                                global=T)

split_list <- split(fttResults1, fttResults1$P3)
set1 <- as.data.frame(split_list[[1]])
set2 <- as.data.frame(split_list[[2]])

# set1 = P1, P2, Others, P3
# set2 = P1, P2, P3, Others

set1$ABBA   # 237
set1$BABA   # 230

set2$ABBA  # 730
set2$BABA  # 230

set1$X2_P # 2.02e-01
set2$X2_P # 1.50e-46

set1$D_jEstimate  # 0.015
set2$D_jEstimate  # 0.521

set1$Fd_1DD4_jEstimate  # 0.002
set1$Fd_D2D4_jEstimate  # 0

set2$Fd_1DD4_jEstimate  # 0.161
set2$Fd_D2D4_jEstimate  # 0

#############################
### Maui Nui 1
Analysis1 <- HC$new("Maui-Nui-concat.NOGAPS.fa")

populations1 <- list(
  P1 = c('Cgyi2'),
  P2 = c('Cgna2','Cmacr'),
  P3 = c('CmXg'),
  Others = c('ChXg','Chaw2','Chaw1')
)

Analysis1$setPopulations(populations1)

Analysis1$prepareFourTaxonTests()

popCombos <- list(
  c(P1 = "P3", P2 = "P2", P3 = "P1", A = "Others")
)

Analysis1$prepareFourTaxonTests(popCombos)

Analysis1$runFourTaxonTests(selections = 'ALL', 
                            blockLength = 1000L)

fttResults1 <- Analysis1$tabulateFourTaxonTests(selections='TESTED',
                                                neat=T,
                                                global=T)

split_list <- split(fttResults1, fttResults1$P1)
set1 <- as.data.frame(split_list[[1]])
set2 <- as.data.frame(split_list[[2]])

# set1 = CmXg, Group1, Group2, Others
# set2 = Group2, Group1, CmXg, Others

set1$ABBA   # 256
set1$BABA   # 295

set2$ABBA  # 307
set2$BABA  # 295

set1$X2_P # 4.587686e-76
set2$X2_P # 2.100774e-47

set1$D_jEstimate  # -0.0698941
set2$D_jEstimate  # 0.02114937

set1$Fd_1DD4_jEstimate  # 0
set1$Fd_D2D4_jEstimate  # 0.0168539

set2$Fd_1DD4_jEstimate  # 0.002355551
set2$Fd_D2D4_jEstimate  # 0
