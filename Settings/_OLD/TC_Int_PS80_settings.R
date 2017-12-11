source(file="Carriers/Simplified/General_settings.R")

#Ac225 recoil, 80 nm PS
#       {Ac,  Fr,  At,  Bi,  Po,  Tl,  Pb,  Bi-stable}
R    = c(0,   0,   0,   0,   0,   0,   0,   0, #Ac
         .65, 0,   0,   0,   0,   0,   0,   0, #Fr
         0,   0,   0,   0,   0,   0,   0,   0, #At
         0,   0,   .43, 0,   0,   0,   0,   0, #Bi
         0,   0,   0,   0,   0,   0,   0,   0, #Po
         0,   0,   0, .45,   0,   0,   0,   0, #Tl
         0,   0,   0,   0, .45,   0,   0,   0, #Pb
         0,   0,   0,   0,   0,   0,   0,   0) #Bi-stable
R    = matrix(R, ncol = length(tau), byrow = T)

# {B, LI, KI, SP, HE, LU, PA, IN, ST, TE, ME, BO, BM, BR, TU, WB, CL}
N[,,] = 0

#Number of carriers over time
PS[,15] = A0*tau[1]/log(2) #tumour

#Number of isotopes per carrier over time
PSi[1,] = c(1, 0, 0, 0, 0, 0, 0, 0)
PSi[1,] = PSi[1,]/sum(PSi[1,])

res  = Carrier.concentration(PSi, NR, tau, dfr, R, t_seq)
PSi  = res$PSi
NR   = res$NR

# plot(t_seq/3600, PS[,1], type="l", xlab="h", ylab="Normalized carrier content")
# plot(t_seq/3600, rowSums(NR), type="l", xlab="h", ylab="Normalized recoil production/s")

#Carrier biostability (half-life)
PSdeg= Inf

#Carrier internalization
PSint= 1