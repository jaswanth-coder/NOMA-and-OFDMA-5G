import numpy as np
import numpy.random as nr
import matplotlib.pyplot as plt
#varinace for user 1 and user 2
delta1sq=1; delta2sq=5;
#delta3sq=(delta1sq*delta2sq)/(delta1sq+delta2sq);     
#creating a SNRdb and converting it into linear scale
SNRdB=np.arange(0,21); SNR=10**(SNRdB/10);    
#a1 power assocaited with user 1
a1=0.9; a2=1-a1;
#it refers to the blocklenghth
blklen=1000000;
#target rate for user1 and user 2
tildeR_1=1; tildeR_2=1; 
#R1 is the thershold rate for target rate snr
R1=2**tildeR_1-1;
R2=2**tildeR_2-1; 
phi=max(R2/a2,R1/(a1-a2*R1));#this will be used in theortical calulation

#Probality of outage we are initalizing
Pout1 = np.zeros(len(SNRdB));
Pout2 = np.zeros(len(SNRdB));
#initalizing the Probability theoroticaly
Pout1_theory = np.zeros(len(SNRdB));
Pout2_theory = np.zeros(len(SNRdB));

#Getting the values and storing the values in the array
for ix in range(len(SNRdB)):
    print(ix);
    
    rhos=SNR[ix];
    #channel fading co-efficents
    h1=np.sqrt(delta1sq/2)*(nr.normal(0,1,blklen)+1j*nr.normal(0,1,blklen))
    h2=np.sqrt(delta2sq/2)*(nr.normal(0,1,blklen)+1j*nr.normal(0,1,blklen))
    #instaneous channel gain ...
    beta_1=np.absolute(h1)**2
    beta_2=np.absolute(h2)**2
    gamma1_u1=(a1*rhos*beta_1)/(a2*rhos*beta_1+1)
    gamma1_u2=(a1*rhos*beta_1)/(a2*rhos*beta_2+1)    
    gamma2_u2=(a2*rhos*beta_2)
    #probalilty
    Pout1[ix]=np.sum(np.log2(1+gamma1_u1)<tildeR_1)/blklen
    Pout2[ix]=np.sum(np.logical_or(np.log2(1+gamma1_u2)<tildeR_1,np.log2(1+gamma2_u2)<tildeR_2))/blklen
    #probability 
    Pout1_theory[ix]=1-np.exp(-R1/(delta1sq*rhos*(a1-a2*R1)));
    Pout2_theory[ix]=1-np.exp(-phi/(rhos*delta2sq));

        
  
#we are ploting the curves
plt.yscale('log')
plt.plot(SNRdB, Pout1,'r-');
plt.plot(SNRdB, Pout1_theory,'rs');
plt.plot(SNRdB, Pout2,'g-.');
plt.plot(SNRdB, Pout2_theory,'go');
plt.grid(1,which='both')
plt.legend(["Outage 1 (Sim)", "Outage 1 (Theory)","Outage 2 (Sim)", "Outage 2 (Theory)"], loc ="lower left");
plt.suptitle('Pout vs SNR for fixed NOMA')
plt.ylabel('Probability of Outage')
plt.xlabel('SNR (dB)') 

