library(pluseTD)
TimeGrid <- c(0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180)
tL <- 10
load(file.path(system.file(package="pluseTD"),'data','rpkmTRUE.RData'))
labexon = rpkmTRUE$labexon
labintr = rpkmTRUE$labintr
totexon = rpkmTRUE$totexon
totintr = rpkmTRUE$totintr
##
load(file.path(system.file(package="pluseTD"),'data','rpkmSim.RData'))

sim_TL_sample = rpkmSim$labexon
sim_PT_sample = rpkmSim$totintr
sim_TT_sample = rpkmSim$totexon

########
#poly
reg2 = regress(sim_TL_sample, sim_TT_sample, sim_PT_sample, TimeGrid, tL, n=2)
reg3 = regress(sim_TL_sample, sim_TT_sample, sim_PT_sample, TimeGrid, tL, n=3)
reg4 = regress(sim_TL_sample, sim_TT_sample, sim_PT_sample, TimeGrid, tL, n=4)
save(reg2, file = "reg2.RData")
save(reg3, file = "reg3.RData")
save(reg4, file = "reg4.RData")
#solver
library(INSPEcT)
mycerIds1 <- newINSPEcT(TimeGrid, tL, labexon[rownames(sim_TL_sample),], totexon[rownames(sim_TL_sample),],
                        labintr[rownames(sim_TL_sample),], totintr[rownames(sim_TL_sample),], BPPARAM=SerialParam())
myat = as.matrix(ratesFirstGuess(mycerIds1, 'synthesis'))
mybt = as.matrix(ratesFirstGuess(mycerIds1, 'degradation'))
myct = as.matrix(ratesFirstGuess(mycerIds1, 'processing'))
solver = list(a = myat, b=mybt, c=myct)
save(solver, file = "solver.RData")

######################
load(file.path(system.file(package="pluseTD"),'data','pluse5.RData'))
load(file.path(system.file(package="pluseTD"),'data','pluse7.RData'))
load(file.path(system.file(package="pluseTD"),'data','pluse9.RData'))
load(file.path(system.file(package="pluseTD"),'data','pluse11.RData'))
##
genename = pluse5@genenames
####################
error_5 = errorMat(pluse5, sim_TL_sample, sim_TT_sample, sim_PT_sample, genename, TimeGrid, tL, idx=1:5, idx2 = 6:13)
error_7 = errorMat(pluse7, sim_TL_sample, sim_TT_sample, sim_PT_sample, genename, TimeGrid, tL, idx=1:7, idx2 = 8:13)
error_9 = errorMat(pluse9, sim_TL_sample, sim_TT_sample, sim_PT_sample, genename, TimeGrid, tL, idx=1:9, idx2 = 10:13)
error_11= errorMat(pluse11, sim_TL_sample, sim_TT_sample, sim_PT_sample, genename, TimeGrid, tL, idx=1:11, idx2 = 12:13)
save(error_5, file = "error_5.RData")
save(error_7, file = "error_7.RData")
save(error_9, file = "error_9.RData")
save(error_11, file = "error_11.RData")
#####rates error###########
load(file.path(system.file(package="pluseTD"),'data','ratesTRUE.RData'))
resmodel20 = correctionParams(resmodel20)
#genename = resmodel20@genenames
same_a = getRates(resmodel20, 'transcription', timevector = TimeGrid)[genename,]
same_b = getRates(resmodel20, 'degradation', timevector = TimeGrid)[genename,]/as.matrix(totexon[genename,]-totintr[genename,])
same_c = getRates(resmodel20, 'processing', timevector = TimeGrid)[genename,]/as.matrix(totintr[genename,])
##
genename=pluse5@genenames[-c(pluse5@fitfailure)]
rates_error_5 = RatesErrorMat(pluse5, sim_TL_sample, sim_TT_sample, sim_PT_sample,
                              same_a,same_b,same_c,
                              genename, TimeGrid, tL, idx=1:5, idx2 = 6:13)
save(rates_error_5, file = "rates_error_5.RData")
##
genename=pluse7@genenames[-c(pluse7@fitfailure)]
rates_error_7 = RatesErrorMat(pluse7, sim_TL_sample, sim_TT_sample, sim_PT_sample,
                              same_a,same_b,same_c,
                              genename, TimeGrid, tL, idx=1:7, idx2 = 8:13)
save(rates_error_7, file = "rates_error_7.RData")
##
genename=pluse9@genenames[-c(pluse9@fitfailure)]
rates_error_9 = RatesErrorMat(pluse9, sim_TL_sample, sim_TT_sample, sim_PT_sample,
                              same_a,same_b,same_c,
                              genename, TimeGrid, tL, idx=1:9, idx2 = 10:13)
save(rates_error_9, file = "rates_error_9.RData")
##
genename=pluse11@genenames[-c(pluse11@fitfailure)]
rates_error_11=RatesErrorMat(pluse11, sim_TL_sample, sim_TT_sample, sim_PT_sample,
                             same_a,same_b,same_c,
                             genename, TimeGrid, tL, idx=1:11, idx2 = 12:13)
save(rates_error_11, file = "rates_error_11.RData")


