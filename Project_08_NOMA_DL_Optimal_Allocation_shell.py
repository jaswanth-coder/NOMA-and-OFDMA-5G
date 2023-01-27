import numpy as np
import numpy.random as nr
import matplotlib.pyplot as plt

#find the optimal power just that sum rate for both user is min
#Avg channel gain
delta1sq=1; delta2sq=5;
#range of SNR and converting it into linear
SNRdB=np.arange(10,31,2); SNR=10**(SNRdB/10);

blklen=1000000;
tildeR_1=1; tildeR_2=1;
#are thershold for u1 and u2
R1=2**(tildeR_1)-1; 
R2=2**(tildeR_2)-1;
#Suboptimal intizational
sr_opt = np.zeros(len(SNRdB));
sr_sub = np.zeros(len(SNRdB));

for ix in range(len(SNRdB)):
    print(ix);
    rhos=SNR[ix];
    #are inizaliing for channel gaussian 
    h1=np.sqrt(delta1sq/2)*(nr.normal(0,1,blklen)+1j*nr.normal(0,1,blklen));
    h2=np.sqrt(delta2sq/2)*(nr.normal(0,1,blklen)+1j*nr.normal(0,1,blklen)); 
    
    #min of their gain is bet1 and maxi ofgain of u2
    beta_1 = np.minimum(np.absolute(h1)**2,np.absolute(h2)**2);
    beta_2 = np.maximum(np.absolute(h1)**2,np.absolute(h2)**2);
    
    #To satisfy QoS constraint at two users
    #
    a2_max=(rhos*beta_1-R1)/(rhos*beta_1*(1+R1));
    a2_min=R2/(rhos*beta_2)
    
    #it is a logical operation. it is the feasible region for optimal allocation
    feasible=a2_max>=a2_min
    
    
    a2_opt=a2_max*feasible 
    a1_opt=(1-a2_opt)*feasible
    
    gamma1_u1_opt=(a1_opt*rhos*beta_1)/(a2_opt*rhos*beta_1+1)
    gamma2_u2_opt=a2_opt*rhos*beta_2
    
    sr_opt[ix]=np.sum(np.log2(1+gamma1_u1_opt)+np.log2(1+gamma2_u2_opt))/np.sum(feasible)
    
    a2=(a2_max+a2_min)/2*feasible
    a1=(1-a2)*feasible 
    
    gamma1_u1_subopt=(a1*rhos*beta_1)/(a2*rhos*beta_1+1)
    gamma2_u2_subopt=a2*rhos*beta_2
    
    sr_sub[ix]=np.sum(np.log2(1+gamma1_u1_subopt)+np.log2(1+gamma2_u2_subopt))/np.sum(feasible)
    
plt.plot(SNRdB, sr_sub,'g-s');
plt.plot(SNRdB, sr_opt,'bo-.');
plt.grid(1,which='both')
plt.legend(["Suboptimal", "Optimal"], loc ="upper left");
plt.suptitle('Sum-Rate with Optimal Power Allocation')
plt.ylabel('Sum-Rate')
plt.xlabel('SNRdB') 

'''observation : suboptimal and optimal line will be nearly same at high snr and at low snr they will have some diffrenec'''
