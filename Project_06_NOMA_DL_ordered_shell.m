
clear all
close all

delta1sq=1; delta2sq=5;
delta3sq = delta1sq*delta2sq/(delta1sq+delta2sq);
SNRdB=0:20; SNR=10.^(SNRdB/10); % SNR range

a1=0.9; a2=0.1;

blklen=1000000; % Block length
tildeR_1=1; tildeR_2=1; % Desired data rates
R1=2^(tildeR_1)-1;
R2=2^(tildeR_2)-1;
phi=max((R2/a2),(R1/(a1-a2*R1)));

for ix=1:length(SNR)
    ix
    outage1=0; outage2=0;5
    
    rhos=SNR(ix);
    h1=sqrt(delta1sq/2)*(randn(1,blklen)+1j*randn(1,blklen));
    h2=sqrt(delta2sq/2)*(randn(1,blklen)+1j*randn(1,blklen));
    beta_1=min(abs(h1).^2,abs(h2).^2);
    beta_2=max(abs(h1).^2,abs(h2).^2);

    gamma1_u1=(a1*rhos*beta_1)./(a2*rhos*beta_1+1);
    gamma1_u2=(a1*rhos*beta_2)./(a2*rhos*beta_2+1);
    gamma2_u2=a2*rhos*beta_2;

    outage1=(log2(1+gamma1_u1)<tildeR_1);%For simulated values
    outage2=((log2(1+gamma1_u2)<tildeR_1)|(log2(1+gamma2_u2)<tildeR_2));

    Pout1(ix)=sum(outage1)/blklen;
    Pout2(ix)=sum(outage2)/blklen;
    
end

Pout1_theory = 1-exp(-R1./(delta3sq*SNR*(a1-a2*R1)));
Pout2_theory = 1-exp(-phi./(SNR*delta1sq))-exp(-phi./(SNR*delta2sq))+exp(-phi./(SNR*delta3sq));


semilogy(SNRdB,Pout1,'r -','LineWidth',2.0)
hold on
semilogy(SNRdB,Pout2,'m -','LineWidth',2.0)
semilogy(SNRdB,Pout1_theory,'r s','LineWidth',2.0,'markerfacecolor','r')
semilogy(SNRdB,Pout2_theory,'m o','LineWidth',2.0,'markerfacecolor','m')
grid on
legend('Outage User 1 (Sim.)','Outage User 2 (Sim.)','Outage User 1 (Theory)','Outage User 2 (Theory)')
xlabel('SNR (dB)')
ylabel('Probability of Outage')
title('Pout vs SNR for ordered NOMA')
