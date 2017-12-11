source(file="Carriers/Simplified/General_settings.R")

#Ac225 recoil loss, 80 nm PS
#       {Ac,  Fr,  At,  Bi,  Po,  Tl,  Pb,  Bi-stable}
R    = c(0,   0,   0,   0,   0,   0,   0,   0, #Ac
         .65, 0,   0,   0,   0,   0,   0,   0, #Fr
         0,   0,   0,   0,   0,   0,   0,   0, #At
         0,   0,   1,   0,   0,   0,   0,   0, #Bi
         0,   0,   0,   0,   0,   0,   0,   0, #Po
         0,   0,   0,   0,   0,   0,   0,   0, #Tl
         0,   0,   0,   0,   0,   0,   0,   0, #Pb
         0,   0,   0,   0,   0,   0,   0,   0) #Bi-stable
R    = matrix(R, ncol = length(tau), byrow = T)

N[,,] = 0

#      {B,  LI,    KI,   SP,    HE, LU,   PA,   IN,   ST, TE, ME,    BO,   BM,    BR, TU, WB, CL}
PSB = c(0,  17.44, 0.99, 17.97, 0,  6.04, 0.15, 1.50, 0,  0,  42.47, 0.43, 13.04, 0,  0,  0,  0)
PSB = PSB/sum(PSB)

#Number of carriers over time
PS[,1]  = A0*tau[1]/log(2)*exp(-0.005/60*t_seq) #Blood
PS[,2]  = (A0*tau[1]/log(2) - PS[,1]) * PSB[2]  #Liver
PS[,3]  = (A0*tau[1]/log(2) - PS[,1]) * PSB[3]  #Kidney
PS[,4]  = (A0*tau[1]/log(2) - PS[,1]) * PSB[4]  #Spleen
PS[,5]  = (A0*tau[1]/log(2) - PS[,1]) * PSB[5]  #Heart
PS[,6]  = (A0*tau[1]/log(2) - PS[,1]) * PSB[6]  #Lungs
PS[,7]  = (A0*tau[1]/log(2) - PS[,1]) * PSB[7]  #Pancreas
PS[,8]  = (A0*tau[1]/log(2) - PS[,1]) * PSB[8]  #Intestines
PS[,9]  = (A0*tau[1]/log(2) - PS[,1]) * PSB[9]  #Stomach
PS[,10] = (A0*tau[1]/log(2) - PS[,1]) * PSB[10] #Testes
PS[,11] = (A0*tau[1]/log(2) - PS[,1]) * PSB[11] #Muscle
PS[,12] = (A0*tau[1]/log(2) - PS[,1]) * PSB[12] #Bone
PS[,13] = (A0*tau[1]/log(2) - PS[,1]) * PSB[13] #Bone Marrow
PS[,14] = (A0*tau[1]/log(2) - PS[,1]) * PSB[14] #Brain
PS[,15] = (A0*tau[1]/log(2) - PS[,1]) * PSB[15] #Tumour
PS[,16] = (A0*tau[1]/log(2) - PS[,1]) * PSB[16] #Whole body
PS[,17] = (A0*tau[1]/log(2) - PS[,1]) * PSB[17] #Clearance


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