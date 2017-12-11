Carrier.concentration <- function(PSi, NR, tau, dfr, R, t_seq) {
  for (t in 2:length(t_seq)) {
    dt = t_seq[t] - t_seq[t-1]
    
    PSi[t,] = PSi[t-1,]*2^(-dt/tau)
    NR[t,]  = NR[t,]  + colSums(PSi[t-1,]*(1 - 2^(-dt/tau))*t(dfr[,])*t(R[,]))
    PSi[t,] = PSi[t,] + colSums(PSi[t-1,]*(1 - 2^(-dt/tau))*t(dfr[,])*t((1 - R[,])))
  }
  
  return(list("PSi"=PSi, "NR"=NR))
}

require(ggplot2)
require(Cairo)
require(progress)

t_seq   = seq(0, tmax, dt)
clearph = 3600*3 #every 3 hours

#Isotope Names
IsoNames = c("Ac225", "Fr221", "At217", "Bi213", "Po213", "Tl209", "Pb209", "Bi-Stable")

#Ac225 half-lifes
tau  = c(10.0*24*3600,
         4.8*60,
         32.3E-6,
         45.6*60,
         4.2E-6,
         2.2*60,
         3.25*3600,
         Inf)

#Ac225 mean energies
Engy = c(5.94,
         6.46,
         7.20,
         1.422*.98 + 5.98*0.02,
         8.53,
         3.98,
         .64,
         0)

#Ac225 chain
#{Ac, Fr, At, Bi, Po, Tl, Pb, Pb-stable}
dfr  = c(0,  0,  0,  0,   0,   0,   0,  0,
         1,  0,  0,  0,   0,   0,   0,  0,
         0,  1,  0,  0,   0,   0,   0,  0,
         0,  0,  1,  0,   0,   0,   0,  0,
         0,  0,  0,  .98, 0,   0,   0,  0,
         0,  0,  0,  0,   .02, 0,   0,  0,
         0,  0,  0,  0,   1,   1,   0,  0,
         0,  0,  0,  0,   0,   0,   1,  0)
dfr  = matrix(dfr, ncol = length(tau), byrow = T)

#Internal decay probabilities
ICP  = c(.00,   #Ac
         .73, #Fr
         1,   #At
         .85, #Bi
         1,   #Po
         1,   #Tl
         .90, #Pb
         1)   #Pb-stable

#Retention due to fast decay probabilities
RDC  = c(0, #Ac
         0, #Fr
         1, #At
         0, #Bi
         1, #Po
         0, #Tl
         0, #Pb
         1) #Pb-stable

#Isotope Colors
IsoCols = c("Red", "Yellow", "Darkblue", "Lightblue", "Gray", "Orange", "Green")

#Organ Names
OrgNames = c("Adrenal","Blood","Bone","Brain","Heart","Intestine","Kidney","Liver","Lung","Muscle","Pancreas",
             "Prostate","RBM","Salivatory","Spleen","Stomach","Thyroid","Tumour","Whole Body","Clearance","Cleared")

Id_WholeBody = match("Whole Body",OrgNames)
Id_Clearance = match("Clearance",OrgNames)
Id_Cleared   = match("Cleared",OrgNames)
Id_Tumour    = match("Tumour",OrgNames)
Id_Blood     = match("Blood",OrgNames)
Id_Bone      = match("Bone",OrgNames)
Id_RBM       = match("RBM",OrgNames)

#            Ad, Bl,   Bo,   Br,   He,  In,   Ki,  Li,   Lu,   Mu,    Pa,  Pr, RBM,  Sa, Sp,  St,  Th, Tu, WB, Cl,  Ce
OrgM     = c(14, 5500, 5000, 1400, 330, 1000, 310, 1800, 1000, 30000, 100, 13, 1500, 85, 180, 150, 20, 10, 0,  300, 300)
OrgM[Id_WholeBody] = 73000 - sum(OrgM[1:Id_WholeBody])
OrgV     = OrgM*1
N        = array(0, dim = c(length(t_seq), length(OrgM), length(tau)))
A        = array(0, dim = c(length(t_seq), length(OrgM), length(tau))) #Activity

RMBLR    = 0.36

#Biodistribution
#           Ad,   Bl,     Bo,    Br,   He,   In,   Ki,   Li,    Lu,   Mu,    Pa,   Pr,   RBM,  Sa,   Sp,   St,   Th,   Tu,   WB,    Cl,    Ce
BioD    = c(0.00, 0.08,   11.52, 0.00, 0.00, 0.00, 2.37, 37.30, 0.00,  0.00, 0.00, 0.00, 0.00, 0.00, 0.17, 0.00, 0.00, 0.00, 37.21, 11.20, 0.00, #Ac
            0.00, 0.92,   5.10,  0.02, 0.65, 6.42, 1.63,  5.41, 0.38,  8.29, 0.16, 0.00, 0.00, 0.00, 0.21, 0.31, 0.00, 0.00, 32.65, 37.74, 0.00, #Fr
            0.00, 100.00, 0.00,  0.00, 0.00, 0.00, 0.00,  0.00, 0.00,  0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,  0.00,  0.00, 0.00, #At
            0.00, 1.92,   15.29, 0.01, 0.02, 0.63, 7.05,  2.12, 0.17,  1.60, 0.03, 0.00, 0.00, 0.00, 0.06, 0.08, 0.00, 0.00, 32.74, 38.24, 0.00, #Bi
            0.00, 100.00, 0.00,  0.00, 0.00, 0.00, 0.00,  0.00, 0.00,  0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,  0.00,  0.00, 0.00, #Po
            0.00, 1.08,   8.40,  0.04, 0.25, 3.71, 2.05,  2.81, 0.26, 10.11, 0.17, 0.00, 0.00, 0.00, 0.21, 0.26, 0.00, 0.00, 22.86, 47.70, 0.00, #Tl
            0.00, 10.89,  3.83,  0.01, 0.06, 1.64, 1.63, 15.96, 0.27,  1.67, 0.00, 0.00, 0.00, 0.00, 0.24, 0.08, 0.00, 0.00, 25.90, 37.79, 0.00, #Pb
            0.00, 100.00, 0.00,  0.00, 0.00, 0.00, 0.00,  0.00, 0.00,  0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,  0.00,  0.00, 0.00) #Pb-stable
BioD    = matrix(BioD, nrow = length(tau), byrow = T)
BioD[,] = BioD[,]/rowSums(BioD[,])

#Carrier internalization
PSint  = array(0, dim = c(length(OrgM)))

#Number of carriers over time
PS   = array(0, dim = c(length(t_seq), length(OrgM)))

#Number of isotopes per carrier over time
PSi  = array(0, dim = c(length(t_seq), length(tau)))

#Number of recoils produced over time
NR   = array(0, dim = c(length(t_seq), length(tau)))

#Storage for decays
DecPS = array(0, dim = c(length(t_seq), length(OrgM), length(tau))) #Carrier decays
DecN  = array(0, dim = c(length(t_seq), length(OrgM), length(tau))) #Free isotope decays

#Decays -> MBq-hr
DMBqHr = 1/(1E6 * 3600)

MBqHrMbq = array(0, dim = c(length(OrgM), length(tau)))

#Dose
DosePS = array(0, dim = c(length(OrgM), length(tau)))
DoseN  = array(0, dim = c(length(OrgM), length(tau)))
MeVtoJoule = 1.60218e-13

#Progress bar
pb = progress_bar$new(total = length(t_seq) - 1)