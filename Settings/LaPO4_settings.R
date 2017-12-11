source(file="Carriers/Simplified/General_settings.r")

#Ac225 recoil production
#       {Ac,  Fr,  At,  Bi,  Po,  Tl,  Pb,  Bi-stable} -> producer
R    = c(0,   0,   0,   0,   0,   0,   0,   0, #Ac
         .5,  0,   0,   0,   0,   0,   0,   0, #Fr
         0,   .5,  0,   0,   0,   0,   0,   0, #At
         0,   0,   .5,  0,   0,   0,   0,   0, #Bi
         0,   0,   0,   0,   0,   0,   0,   0, #Po
         0,   0,   0,   .5,  0,   0,   0,   0, #Tl
         0,   0,   0,   0,   .5,  0,   0,   0, #Pb
         0,   0,   0,   0,   0,   0,   0,   0) #Bi-stable
#50% retention of daughters (see paper)
R    = matrix(R, ncol = length(tau), byrow = T)

N[,,] = 0

#Number of carriers over time
PS[,Id_Blood]                     =  A0*tau[1]/log(2) * 2^(-t_seq / (8*3600)) #Blood
PS[,match("Liver",OrgNames)]      = (A0*tau[1]/log(2) - PS[,Id_Blood]) * 34/100
PS[,match("Kidney",OrgNames)]     = (A0*tau[1]/log(2) - PS[,Id_Blood]) * 8/100
PS[,match("Spleen",OrgNames)]     = (A0*tau[1]/log(2) - PS[,Id_Blood]) * 1/100
PS[,match("Heart",OrgNames)]      = (A0*tau[1]/log(2) - PS[,Id_Blood]) * 2/100
PS[,match("Lung",OrgNames)]       = (A0*tau[1]/log(2) - PS[,Id_Blood]) * 5/100
PS[,match("Intestine",OrgNames)]  = (A0*tau[1]/log(2) - PS[,Id_Blood]) * 13/100
PS[,match("Stomach",OrgNames)]    = (A0*tau[1]/log(2) - PS[,Id_Blood]) * 0/100
PS[,match("Muscle",OrgNames)]     = (A0*tau[1]/log(2) - PS[,Id_Blood]) * 15/100
PS[,match("Whole Body",OrgNames)] = (A0*tau[1]/log(2) - PS[,Id_Blood]) * 22/100


#Number of isotopes per carrier over time
PSi[1,] = c(1, 0, 0, 0, 0, 0, 0, 0)
PSi[1,] = PSi[1,]/sum(PSi[1,])

res  = Carrier.concentration(PSi, NR, tau, dfr, R, t_seq)
PSi  = res$PSi
NR   = res$NR

#Carrier biostability (half-life)
PSdeg= Inf