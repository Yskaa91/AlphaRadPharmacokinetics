if (!exists("cl")) {
  cat("\014")
  rm(list = ls())
}

ptm <- proc.time()

Tprefix = "7d"

# time interval
dt = 60

# total time to run
if (Tprefix == "4h") {
  tmax = 3600*4
} else if (Tprefix == "24h") {
  tmax = 3600*24
} else {
  Tprefix = "7d"
  tmax = 3600*24*7
}

NonDecayCorrection = F
  
# Initial activity in Bq (assumes 1 isotope per carrier)
A0   = 1E6 # 1 MBq

plot_min = 1E-10
plot_max = 1E5

doplot = F

# for (prefix in c("HUM195", "Liposomes", "Polymersomes", "Trastuzumab", "LaPO4", "LaPO4_Gold", "Theoretical", "Theoretical2", "PSMA")) {
for (prefix in c("TheoreticalPSMA", "TheoreticalTRAS", "TheoreticalPLMS", "TheoreticalLPS", "TheoreticalLAPO4", "TheoreticalGOLD")) {
# for (prefix in c("Liposomes", "Polymersomes")) {
  source(file=paste("Simplified/Settings/", prefix, "_settings.R", sep=""), local = TRUE)
  
  ## Calculate free isotopes
  caprod = 0
  for (t in 2:length(t_seq)) {
    isopr = array(0, dim = length(tau))
    for (k in 1:length(OrgM)) {
        N[t,k,] = N[t-1,k,]*2^(-(t_seq[t] - t_seq[t-1])/tau)           #Decay of the isotope
        A[t,k,] = (N[t-1,k,] - N[t,k,])/(t_seq[t] - t_seq[t-1])
        
        caprod   = NR[t-1,]*PS[t-1,k] +                                #Recoil production
          PSi[t-1,]*PS[t-1,k]*(1 - 2^(-(t_seq[t] - t_seq[t-1])/PSdeg)) #Carrier decay
        
        fiprod = colSums(N[t-1,k,]*                                    #Isotope production from mother
                    (1 - 2^(-(t_seq[t] - t_seq[t-1])/tau))*t(dfr))
        
        # browse
        N[t,k,] = N[t,k,] + (fiprod + caprod)*RDC                      #Fast decay
	      fiprod  = fiprod*(1 - RDC)
	      caprod  = caprod*(1 - RDC)

        if (k == Id_Blood || k == Id_Bone)
          isopr = isopr + fiprod + caprod                              #Full release
        else if (k == Id_Clearance)
          N[t,k,] = N[t,k,] + fiprod + caprod                          #Full retention
        else {
          #PSint -> Fraction of carrier internalized
          N[t,k,] = N[t,k,] + (fiprod + caprod*PSint[k])*ICP           #Retention probability
          isopr   = isopr + (fiprod + caprod*PSint[k])*(1 - ICP) +     #Release probability
                    caprod*(1 - PSint[k])
        }
    }
    
    for (k in 1:length(OrgM))                                          #Distribute produced free isotopes
      N[t,k,] = N[t,k,] + isopr*BioD[,k]                               #according to biodistribution
    
    if (t_seq[t] %% clearph == 0) {                                    #Regular clearance flush
      N[t,Id_Cleared,]   = N[t,Id_Cleared,] + N[t,Id_Clearance,]
      N[t,Id_Clearance,] = 0
    }
    
    pb$tick()
  }
  
  # tauM <- matrix(rep(tau, length(OrgM)), ncol=length(tau), byrow = T)
  
  ## Calculate decays
  for (t in 2:length(t_seq))
    for (k in 1:length(OrgM)) {
      DecPS[t,k,] = PSi[t-1,]*PS[t-1,k]*(1 - 2^(-(t_seq[t] - t_seq[t-1])/tau))
      #[decays in carrier]
      DecN[t,k,]  = N[t-1,k,]*(1 - 2^(-(t_seq[t] - t_seq[t-1])/tau))
      #[decays in free isotopes]
    }
  
  #Transfer decays of 1500 ml of blood to bone marrow
  DecPS[,Id_RBM,]    = DecPS[,Id_RBM,] + DecPS[,Id_Blood,]* RMBLR*OrgV[Id_RBM] /OrgV[Id_Blood]
  DecN [,Id_RBM,]    = DecN[,Id_RBM,]  + DecN[,Id_Blood,] * RMBLR*OrgV[Id_RBM] /OrgV[Id_Blood]
  DecPS[,Id_Blood,]  = DecPS[,Id_Blood,]*(OrgV[Id_Blood] -  RMBLR*OrgV[Id_RBM])/OrgV[Id_Blood]
  DecN [,Id_Blood,]  = DecN[,Id_Blood,] *(OrgV[Id_Blood] -  RMBLR*OrgV[Id_RBM])/OrgV[Id_Blood]
  
  ## MBq-hr/MBq
  MBqHrMbq[,] = colSums(DecPS[,,] + DecN[,,])*DMBqHr / (A0/1E6)
  
  ## Dose
  DosePS[,]   = sweep(colSums(DecPS[,,]),MARGIN=2,Engy,`*`)*MeVtoJoule/OrgM*1000 #g to kg
  DoseN[,]    = sweep(colSums(DecN[,,] ),MARGIN=2,Engy,`*`)*MeVtoJoule/OrgM*1000 #g to kg
  
  if (NonDecayCorrection) {
    #Non-decayed particle correction
    DosePS = DosePS/(1 - 2^(-tmax/tau[1]))
    DoseN  = DoseN /(1 - 2^(-tmax/tau[1]))
    
    mbqhmbqcheck = sum(MBqHrMbq[,1]) / ((tau[1] - tau[1]*2^(-tmax/tau[1]))/(3600*log(2))) - 1
    
    normcheckadj = abs(sum(PS[length(t_seq),])*sum(PSi[length(t_seq),]) +  # Residual PS content
                         sum(N[length(t_seq),,]) -                           # Residual free isotopes
                         sum(PS[1,]))/sum(PS[1,])                            # Total injected
    
    cat("Normality check  (", formatC(abs(normcheckadj), digits=3, format="E"), ") : ")
    if (abs(normcheckadj) < 0.01) {
      cat("OK\n")
    } else {
      cat("ERROR\n")
    }
    cat("MBq-hr/MBq check (", formatC(abs(mbqhmbqcheck), digits=3, format="E"), ") : ")
    if (abs(mbqhmbqcheck) < 0.01) {
      cat("OK\n\n")
    } else {
      cat("ERROR\n\n")
    }
  }
  
  options(warn=-1)
  cat("Final isotopes/mL\n")
  if (doplot) {
    for (Norg in 1:(length(OrgM)-1)) {
      miny = N[,Norg,]
      miny = min(miny[miny > 0])
      maxy = max(N[,Norg,])
      if (Norg != Id_Clearance) {
        Cairo(width = 700, height = 400, file=paste('Plots/', Tprefix ,'/', prefix, '_', OrgNames[Norg], '.png', sep=''), type="png")
        plot(t_seq/3600, A[,Norg,1]/OrgV[Norg], log="xy", col=IsoCols[1], ylim=c(plot_min, plot_max), xlim=c(t_seq[2], max(t_seq))/3600,
             type="l", xlab="h", ylab="#kBq/mL", main=OrgNames[Norg], lwd=2)
        grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
        
        for (i in 2:length(IsoCols))
          lines(t_seq/3600, A[,Norg,i]/OrgV[Norg], ylim=c(1E-6, 1E8), xlim=c(t_seq[5], max(t_seq))/3600, col=IsoCols[i], lwd=2)
        
        legend("bottomleft", legend = IsoNames[1:length(tau)-1], col=IsoCols, lty=rep(1,length(tau) - 1), ncol=7, lwd=3)
        dev.off()
        
      } else {
        if (tmax > clearph) {
          c_ind = tail(t_seq %% clearph == 0, -1)
          Cairo(width = 700, height = 400, file=paste('Plots/', Tprefix ,'/', prefix, '_', OrgNames[Norg], '.png', sep=''), type="png")
          plot(t_seq[c_ind]/3600, A[c_ind,Norg,1]/OrgV[Norg], log="xy", col=IsoCols[1], ylim=c(1E-10, 1E5), xlim=c(min(t_seq[c_ind]), max(t_seq))/3600,
               type="p", xlab="h", ylab="#kBq/mL", main=OrgNames[Norg])
          grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
          
          for (i in 2:length(IsoCols))
            lines(t_seq[c_ind]/3600, A[c_ind,Norg,i]/OrgV[Norg], ylim=c(1E-10, 1E5), xlim=c(min(t_seq[c_ind]), max(t_seq))/3600, col=IsoCols[i], type="p")
          
          legend("bottomleft", legend = IsoNames[1:length(tau)-1], col=IsoCols, lty=rep(1,length(tau) - 1), ncol=7)
          dev.off()
        }
      }
      
      cat(OrgNames[Norg] ,": ", formatC(sum(N[length(t_seq),Norg,]/OrgV[Norg]), digits=3, format="E"), "\n")
    }
  }
  
  filename = paste('CSV/', Tprefix ,'_', prefix, '.csv', sep='')
  cat("MBq-Hr/MBq\n\"Organ\",", file = filename, append = FALSE)
  write.table(MBqHrMbq, file = filename, row.names=OrgNames, na="", col.names=IsoNames, sep=",", append="true")
  cat("\n", file = filename, append = TRUE)
  
  cat("MBq-Hr/MBq (%)\n\"Organ\",", file = filename, append = TRUE)
  write.table(rowSums(MBqHrMbq)/sum(MBqHrMbq)*100, file = filename, row.names=OrgNames, na="", col.names=c("Total"), sep=",", append="true")
  cat("\n", file = filename, append = TRUE)
  
  cat("Isotope Doses (Gy)\n\"Organ\",", file = filename, append = TRUE)
  write.table(DosePS+DoseN, file = filename, row.names=OrgNames, na="", col.names=IsoNames, sep=",", append="true")
  cat("\n", file = filename, append = TRUE)
  
  cat("Total Doses (Gy)\n\"Organ\",", file = filename, append = TRUE)
  write.table(t(rbind(rowSums(DosePS),rowSums(DoseN),rowSums(DosePS)+rowSums(DoseN))), file = filename, row.names=OrgNames, na="", col.names=c("Carrier", "Free", "Total"), sep=",", append="true")
}

cat("Time elapsed:", (proc.time() - ptm)[3], "s")
