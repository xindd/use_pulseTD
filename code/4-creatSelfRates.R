library(pluseTD)
# rpkmTRUE
load(file.path(system.file(package="pluseTD"),'data','rpkmTRUE.RData'))
TimeGrid <- c(0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180)
tL <- 10
labexon = rpkmTRUE$labexon
totexon = rpkmTRUE$totexon
totintr = rpkmTRUE$totintr
##
load(file.path(system.file(package="pluseTD"),'data','rpkmSim.RData'))
##
########################
gnum = 1000
sim_TL_sample = rpkmSim$labexon[1:gnum, ]
sim_PT_sample = rpkmSim$totintr[1:gnum, ]
sim_TT_sample = rpkmSim$totexon[1:gnum, ]

#pluseRates = estimateParams(sim_TL_sample, sim_TT_sample, sim_PT_sample, TimeGrid, tL, loopnumber=50)
# = correctionParams(pluseRates)
#save(pluseRates, file = "pluseRates.RData")
load(file.path(system.file(package="pluseTD"),'data','pluseRates.RData'))

simgene = pluseRates@genenames
###########
pluse5 = estimateParams(sim_TL_sample[simgene,1:5],
                        sim_TT_sample[simgene,1:5],
                        sim_PT_sample[simgene,1:5], TimeGrid[1:5], tL, loopnumber=60)
#pluse5 = correctionParams(pluse5)
save(pluse5, file = "pluse5.RData")
pluse7 =  estimateParams(sim_TL_sample[simgene,1:7],
                         sim_TT_sample[simgene,1:7],
                         sim_PT_sample[simgene,1:7], TimeGrid[1:7], tL, loopnumber=60)
#pluse7 = correctionParams(pluse7)
save(pluse7, file = "pluse7.RData")
pluse9 =  estimateParams(sim_TL_sample[simgene,1:9],
                         sim_TT_sample[simgene,1:9],
                         sim_PT_sample[simgene,1:9], TimeGrid[1:9], tL, loopnumber=60)
#pluse9 = correctionParams(pluse9)
save(pluse9, file = "pluse9.RData")
pluse11= estimateParams(sim_TL_sample[simgene,1:11],
                        sim_TT_sample[simgene,1:11],
                        sim_PT_sample[simgene,1:11], TimeGrid[1:11], tL, loopnumber=60)
#pluse11 = correctionParams(pluse11)
save(pluse11, file = "pluse9.RData")
