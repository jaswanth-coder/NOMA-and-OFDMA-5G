import numpy as np
import numpy.random as nr
import matplotlib.pyplot as plt
#SIC is also called the basic b-blast recevier
#in the uplink sys base staion recive will recive both users data
'''' base station deodes the data for The stroger user will decode its data with interference.Power allocatio will be there  '''
#Avg channel gain users are in decreasing order..
delta1sq=5; delta2sq=1;
SNRdB=np.arange(0,31,2); SNR=10**(SNRdB/10);  
#power allocation for u1 and u2  
a1=0.9; a2=1-a1;#This the power controlling mechanisnm to minimize the interference while decoding
 
blklen=1000000;
#traget data rate
tildeR_1=1; tildeR_2=1; 

R1=2**tildeR_1-1;
R2=2**tildeR_2-1;
phi=max(R2/a2,R1/(a1-a2*R1));
#initiazlisation of outage probaility for outage 1 and ioutage 2
Pout1 = np.zeros(len(SNRdB));
Pout2 = np.zeros(len(SNRdB));
Pout1_theory = np.zeros(len(SNRdB));
Pout2_theory = np.zeros(len(SNRdB));

for ix in range(len(SNRdB)):
    print(ix);
    
    rhos=SNR[ix];
    h1=np.sqrt(delta1sq/2)*(nr.normal(0,1,blklen)+1j*nr.normal(0,1,blklen));
    h2=np.sqrt(delta2sq/2)*(nr.normal(0,1,blklen)+1j*nr.normal(0,1,blklen));
    
    #Instaneous channel gain
    beta_1=np.absolute(h1)**2; beta_2=np.absolute(h2)**2;
    #SINR of u1 at bs and SINR of u2 at bs
    #bs frist decode the strong user frist and deocde the weaker user after doing SIC
    gamma1=(a1*rhos*beta_1)/(a2*rhos*beta_2+1)
    gamma2=a2*rhos*beta_2
    
    #outage calculation for u1 and u2
    Pout1[ix]=np.sum(np.log2(1+gamma1)<tildeR_1)/blklen
    #If error has while decoding u1 data thata error can progate into u2 at the base station.
    Pout2[ix]=np.sum(np.logical_or(np.log2(1+gamma1)<tildeR_1,np.log2(1+gamma2)<tildeR_2))/blklen
    

    
    phi2=(R1/(a1*delta1sq))+(R2/(a2*delta2sq))+((R1*R2)/(a1*delta1sq));
    phi3=1+((R1*a2*delta2sq)/(a1*delta1sq));
    #theorotical calculation of the outage.
    Pout1_theory[ix]=1-((np.exp(-R1/(a1*rhos*delta1sq)))/(1+((R1*a2*delta2sq)/(a1*delta1sq))));
    Pout2_theory[ix]=1-((np.exp(-phi2/rhos))/(phi3));


        
  
#Plotting code 
plt.yscale('log')
plt.plot(SNRdB, Pout1,'r-');
plt.plot(SNRdB, Pout1_theory,'rs');
plt.plot(SNRdB, Pout2,'g-.');
plt.plot(SNRdB, Pout2_theory,'go');
plt.grid(1,which='both')
plt.legend(["Outage 1 (Sim)", "Outage 1 (Theory)","Outage 2 (Sim)", "Outage 2 (Theory)"], loc ="upper right");
plt.suptitle('Pout vs SNR for UL NOMA')
plt.ylabel('Probability of Outage')
plt.xlabel('SNR (dB)') 

#observation
'''u2 outage is high and u1 outage is less because error progation maybe occcur, outage will saturate at high SNR'''

