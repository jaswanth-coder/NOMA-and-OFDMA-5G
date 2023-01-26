import numpy as np
import numpy.random as nr
import matplotlib.pyplot as plt

#statistical variance 
delta1sq=1; delta2sq=5;
#delta3 is the harmonic value of delta1 and delta2
delta3sq=(delta1sq*delta2sq)/(delta1sq+delta2sq);     
SNRdB=np.arange(0,21); SNR=10**(SNRdB/10);    
#path co-efficengt for u1 and u2
a1=0.9; a2=1-a1;
 
#blocklength for coding
blklen=1000000;
#traget data rate /user quality data rate.
#in praftical it depend son the what kind of data rate and whta is the user need.
tildeR_1=1; tildeR_2=1; 
#Thershold snr
R1=2**tildeR_1-1;
R2=2**tildeR_2-1;
#phi is used in the theoritacal calculation
phi=max(R2/a2,R1/(a1-a2*R1));

#we are the initalizing Probaliltly of outage over tghe length of SNR 
Pout1 = np.zeros(len(SNRdB));
Pout2 = np.zeros(len(SNRdB));

Pout1_theory = np.zeros(len(SNRdB));
Pout2_theory = np.zeros(len(SNRdB));

for ix in range(len(SNRdB)):
    print(ix);
    
    rhos=SNR[ix];
    #these are the complex channels complex gassuina. mean=0
    h1=np.sqrt(delta1sq/2)*(nr.normal(0,1,blklen)+1j*nr.normal(0,1,blklen))
    h2=np.sqrt(delta2sq/2)*(nr.normal(0,1,blklen)+1j*nr.normal(0,1,blklen))
    
    #we are taking the minimum value and maximum value for u1 and u2
    #This is the diffrenec btw fixed NOMA and ordred NOMA line no(41-42)
    beta_1=np.minimum(np.absolute(h1)**2,np.absolute(h2)**2)
    beta_2=np.maximum(np.absolute(h1)**2,np.absolute(h2)**2)
    
    #The desired signal power for u1,u2,u3 SINR value 
    gamma1_u1=(a1*rhos*beta_1)/(a2*rhos*beta_1+1)
    gamma1_u2=(a1*rhos*beta_2)/(a2*rhos*beta_2+1)
    #after doing the SIC
    gamma2_u2=a2*rhos*beta_2
    
    #we are simulating the outage probability for u1 and u2
    Pout1[ix]=np.sum(np.log2(1+gamma1_u1)<tildeR_1)/blklen
    Pout2[ix]=np.sum(np.logical_or(np.log2(1+gamma1_u2)<tildeR_1,np.log2(1+gamma2_u2)<tildeR_2))/blklen
    
    Pout1_theory[ix]=1-np.exp(-R1/(delta3sq*rhos*(a1-a2*R1)));
    Pout2_theory[ix]=1-np.exp(-phi/(rhos*delta1sq))-np.exp(-phi/(rhos*delta2sq))+ np.exp(-phi/(rhos*delta3sq));

        
  
#ploting the graph section.
plt.yscale('log')
plt.plot(SNRdB, Pout1,'r-');
plt.plot(SNRdB, Pout1_theory,'rs');
plt.plot(SNRdB, Pout2,'g-.');
plt.plot(SNRdB, Pout2_theory,'go');
plt.grid(1,which='both')
plt.legend(["Outage 1 (Sim)", "Outage 1 (Theory)","Outage 2 (Sim)", "Outage 2 (Theory)"], loc ="lower left");
plt.suptitle('Pout vs SNR for Ordered NOMA')
plt.ylabel('Probability of Outage')
plt.xlabel('SNR (dB)') 

#why is there intersection in the graph for sim and theoritical
'''at low snr if decoding of u1 at u2 then error is going to propagate
Then the outage will increse at low snr'''
#this is called SIC failure...

