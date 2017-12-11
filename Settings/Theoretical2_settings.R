source(file="Carriers/Simplified/General_settings.r")

#Ac225 recoil production
#       {Ac,  Fr,  At,  Bi,  Po,  Tl,  Pb,  Bi-stable} -> producer
R    = c(0,   0,   0,   0,   0,   0,   0,   0, #Ac
         1,   0,   0,   0,   0,   0,   0,   0, #Fr
         0,   0,   0,   0,   0,   0,   0,   0, #At
         0,   0,   0,   0,   0,   0,   0,   0, #Bi
         0,   0,   0,   0,   0,   0,   0,   0, #Po
         0,   0,   0,   0,   0,   0,   0,   0, #Tl
         0,   0,   0,   0,   0,   0,   0,   0, #Pb
         0,   0,   0,   0,   0,   0,   0,   0) #Bi-stable
R    = matrix(R, ncol = length(tau), byrow = T)

N[,,] = 0

#Number of carriers over time
PS[,match("Tumour",OrgNames)] = A0*tau[1]/log(2)

#Number of isotopes per carrier over time
PSi[1,] = c(1, 0, 0, 0, 0, 0, 0, 0)
PSi[1,] = PSi[1,]/sum(PSi[1,])

res  = Carrier.concentration(PSi, NR, tau, dfr, R, t_seq)
PSi  = res$PSi
NR   = res$NR

#Carrier biostability (half-life)
PSdeg= Inf

#Carrier internalization
PSint[match("Tumour",OrgNames)] = 0