source(file="Carriers/Simplified/General_settings.r")

#Ac225 recoil production, 80nm polymersomes
#       {Ac,  Fr,  At,  Bi,  Po,  Tl,  Pb,  Bi-stable} -> producer
R    = c(0,   0,   0,   0,   0,   0,   0,   0, #Ac
         .65, 0,   0,   0,   0,   0,   0,   0, #Fr
         0,   0,   0,   0,   0,   0,   0,   0, #At
         0,   0,   .43, 0,   0,   0,   0,   0, #Bi
         0,   0,   0,   0,   0,   0,   0,   0, #Po
         0,   0,   0, .45,   0,   0,   0,   0, #Tl
         0,   0,   0,   0, .45,   0,   0,   0, #Pb
         0,   0,   0,   0,   0,   0,   0,   0) #Bi-stable
R    = matrix(R, ncol = length(tau), byrow = T)

N[,,] = 0

#Number of carriers over time
PS[,Id_Blood]                     =  A0*tau[1]/log(2) * 2^(-t_seq / (16*3600)) #Blood
PS[,match("Tumour",OrgNames)]     = (A0*tau[1]/log(2) - PS[,Id_Blood]) * 100/100


#Number of isotopes per carrier over time
PSi[1,] = c(1, 0, 0, 0, 0, 0, 0, 0)
PSi[1,] = PSi[1,]/sum(PSi[1,])

res  = Carrier.concentration(PSi, NR, tau, dfr, R, t_seq)
PSi  = res$PSi
NR   = res$NR

#Carrier biostability (half-life)
PSdeg= Inf

#Carrier internalization
PSint[match("Tumour",OrgNames)] = 1