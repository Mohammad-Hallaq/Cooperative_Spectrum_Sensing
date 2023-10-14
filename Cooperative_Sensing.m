clc
close all
clear


NSU = 3:30;

pfa_OR = 1 - 0.9.^(1./NSU); 

figure;
plot(NSU,pfa_OR,'-o')
grid on
xlim([3 30])
title("Pfa-OR- Vs n")
xlabel("Number of Secondary Users")
ylabel("Local Probability of False Alarm Pfa")


pfa_AND =  0.1.^(1./NSU); 

figure;
plot(NSU,pfa_AND,'-o')
grid on
xlim([3 30])
title("Pfa-AND- Vs n")
xlabel("Number of Secondary Users")
ylabel("Local Probability of False Alarm Pfa")

pfa_K_N = zeros(1,length(NSU));
for i=1:length(NSU)
    
    min_lim = ceil((NSU(i)+1)/2);
    pfa_K_N(i) = inv_K_N(0.1,NSU(i),min_lim);
    

end

figure;
plot(NSU,pfa_K_N,'-o')
grid on
xlim([3 30])
title("Pfa-Majority- Vs n")
xlabel("Number of Secondary Users")
ylabel("Local Probability of False Alarm Pfa")


%%
%%%%%%% OR RULE EXPERIMENTS
clc
clear


zfun = @(Pfa1,Pfa2) Pfa1 + Pfa2 - Pfa1.*Pfa2 ;
zhandle = fcontour(zfun);

zhandle.LevelList = 0.007;
xlabel Pfa_-ed
ylabel Pfa_-mme
title('Pfa_-ed +(1-Pfa_-ed)*Pfa_-mme')
grid on
zhandle.YRange = [0,0.05];
zhandle.XRange = [0,0.05];
axis equal
%%
clc
clear
close all

Nsu = 3:30;
lim = 1-0.9.^(1./Nsu);
pfa1 = lim;
for i =1:length(Nsu)
    for j =1:length(pfa1)
        
        pfa2(i,j) = 1 - (0.9.^(1/Nsu(i)))/(1-pfa1(j));
        
    end
end

figure;
surf(Nsu,pfa1,pfa2)
xlabel("NSU")
ylabel("PFA1")

%%

clc
clear 
close all

L=8; %Smoothing Factor
M=1; %Number of Antennas
P = M; 
pfa_mme = 0.005;
pfa_ed = (0.007 - pfa_mme)/(1 - pfa_mme); %Probability of False Alarm
%pfa_mme = 0:0.001:0.007; %Probability of False Alarm
N = 1e3:2000:2.2*1e4;%:5000:5*1e4; %Number of Samples
snr_dB = -25:5; % SNR in decibels
snr = 10.^(snr_dB./10); % Linear Value of SNR
%% Simulation to plot Probability of Detection (Pd) vs. Signal to Noise Ratio (SNR) 
%Changing Number of samples

for num=1:length(N)
    
for m = 1:length(snr)
    
    detect = 0;
    thresh_ed = (sqrt(2).*qfuncinv(pfa_ed)./sqrt(N(num)))+ 1; % Theoretical value of threshold
    pfa_mme = (0.007 - pfa_ed)/(1 - pfa_ed);
    part1 = ((sqrt(N(num))+sqrt(M*L))/(sqrt(N(num))-sqrt(M*L)))^2;
    part2 = invtw(1-pfa_mme)*(sqrt(N(num))+sqrt(M*L))^(-2/3);
    part3 = 1 + part2/((N(num)*M*L)^(1/6));
    thresh_mme = part1*part3;  % Theoretical value of threshold
     
for k=1:2000 % Number of Monte Carlo Simulations
    
    n=0:N(num)-1;
    S = 1.*sin(2*pi*1/N(num)*n); %A sinusoidal signal
    Ps = sum(abs(S).^2)/length(S); %Power of the signal
    y = awgn(S,snr_dB(m), 'measured'); %Adding the effect of AWGN

 energy = abs(y).^2; % Energy of received signal over N samples
 Test_ed = (1/N(num)).*sum(energy)/(Ps/(snr(m))); % Test Statistic for the Energy Detection
 
 if(Test_ed >= thresh_ed)  % Check whether the received energy is greater than threshold, if so, increment Pd (Probability of detection) detctor by 1
     
     detect = detect+1;
     
 else %Moving to MME
     
    y_mid = [y zeros(1,L)];
    yy = transpose(y_mid);

    summ = zeros(M*L,M*L);
    for i=1:N(num)+1
        
        y_h =  flip(yy(i:L+i-1,:));
        y_hat = reshape(y_h,[M*L,1]);
        summ = summ + y_hat*y_hat';

    end
    
    R = summ/N(num);  %Statistical Covariance Matrix
    lamda = eig(R); %Eigenvalues of R
    Test_mme = max(lamda)/min(lamda);% Test Statistic for MME
    
    
    if(Test_mme >= thresh_mme)
        
        detect = detect+1;
        
    end
     
 end
 
end
Pd(num,m) = detect/k; %Computing the Probability of Detection
end
end

Pd_tot = 1-(1-Pd).^15;

figure;
plot(snr_dB,Pd_tot,'-o')
grid on
title("Pd Vs Pfa")
xlabel("Probability of False Alarm for Energy Detection Stage")
ylabel("Cooperative Probability of Detction Pd")
legend('N=10000','N=12000','N=14000','N=16000','N=18000','N=20000')

figure;
surf(pfa_ed,N(1:7),Pd_tot(1:7,:))
xlabel("Probability of False Alarm for Energy Detection Stage")
ylabel("Number of Samples")

%%
clc
clear 
close all

L=32; %Smoothing Factor
M=1; %Number of Antennas
P = M; 
pfa_ed = 0:0.005:0.03; %Probability of False Alarm
N = 1000:500:3e3;%:1000:0.9*1e4;%:5000:5*1e4; %Number of Samples
snr_dB = -21; % SNR in decibels
snr = 10^(snr_dB./10); % Linear Value of SNR
%% Simulation to plot Probability of Detection (Pd) vs. Signal to Noise Ratio (SNR) 
%Changing Number of samples

for num=1:length(N)
    
for m = 1:length(pfa_ed)
    
    detect = 0;
for k=1:2000 % Number of Monte Carlo Simulations
    
    n=0:N(num)-1;
    S = 1.*sin(2*pi*1/N(num)*n); %A sinusoidal signal
    Ps = sum(abs(S).^2)/length(S); %Power of the signal
    y = awgn(S,snr_dB, 'measured'); %Adding the effect of AWGN

 energy = abs(y).^2; % Energy of received signal over N samples
 Test_ed = (1/N(num)).*sum(energy)/(Ps/(snr)); % Test Statistic for the Energy Detection
 thresh_ed = (sqrt(2).*qfuncinv(pfa_ed(m))./sqrt(N(num)))+ 1; % Theoretical value of threshold
 
 pfa_mme = (0.0345 - pfa_ed(m))/(1 - pfa_ed(m));
 part1 = ((sqrt(N(num))+sqrt(M*L))/(sqrt(N(num))-sqrt(M*L)))^2;
 part2 = invtw(1-pfa_mme)*(sqrt(N(num))+sqrt(M*L))^(-2/3);
 part3 = 1 + part2/((N(num)*M*L)^(1/6));
 thresh_mme = part1*part3;  % Theoretical value of threshold
 
 if(Test_ed >= thresh_ed)  % Check whether the received energy is greater than threshold, if so, increment Pd (Probability of detection) detctor by 1
     
     detect = detect+1;
     
 else %Moving to MME
     
    y_mid = [y zeros(1,L)];
    yy = transpose(y_mid);

    summ = zeros(M*L,M*L);
    for i=1:N(num)+1
        
        y_h =  flip(yy(i:L+i-1,:));
        y_hat = reshape(y_h,[M*L,1]);
        summ = summ + y_hat*y_hat';

    end
    
    R = summ/N(num);  %Statistical Covariance Matrix
    lamda = eig(R); %Eigenvalues of R
    Test_mme = max(lamda)/min(lamda);% Test Statistic for MME
    
    
    if(Test_mme >= thresh_mme)
        
        detect = detect+1;
        
    end
     
 end
 
end
Pd(num,m) = detect/k; %Computing the Probability of Detection
end
end

Pd_tot = 1-(1-Pd).^3;

figure;
plot(pfa_ed,Pd_tot,'-o')
grid on
title("Pd Vs Pfa")
xlabel("Probability of False Alarm for Energy Detection Stage")
ylabel("Cooperative Probability of Detction Pd")
legend('N=1000','N=1500','N=2000','N=2500','N=3000')

figure;
surf(pfa_ed,N,Pd_tot)
xlabel("Probability of False Alarm for Energy Detection Stage")
ylabel("Number of Samples")
%%
%%%%%%% AND RULE EXPERIMENTS
clc
clear


zfun = @(Pfa1,Pfa2) Pfa1 + Pfa2 - Pfa1.*Pfa2 ;
zhandle = fcontour(zfun);

zhandle.LevelList = 0.4642;
xlabel Pfa_-ed
ylabel Pfa_-mme
title('Pfa_-ed +(1-Pfa_-ed)*Pfa_-mme')
grid on
zhandle.YRange = [0,1];
zhandle.XRange = [0,1];
axis equal
%%
clc
clear 
close all

L=8; %Smoothing Factor
M=1; %Number of Antennas
P = M; 
pfa_ed = 0:0.05:0.45; %Probability of False Alarm
pfa_mme = 0:0.05:0.45; %Probability of False Alarm
N = 2*1e4:2000:4*1e4; %Number of Samples
snr_dB = -21; % SNR in decibels
snr = 10^(snr_dB./10); % Linear Value of SNR
%% Simulation to plot Probability of Detection (Pd) vs. Signal to Noise Ratio (SNR) 
%Changing Number of samples

for num=1:length(N)
    
for m = 1:length(pfa_ed)
    
    detect = 0;
for k=1:2000 % Number of Monte Carlo Simulations
    
    n=0:N(num)-1;
    S = 1.*sin(2*pi*1/N(num)*n); %A sinusoidal signal
    Ps = sum(abs(S).^2)/length(S); %Power of the signal
    y = awgn(S,snr_dB, 'measured'); %Adding the effect of AWGN

 energy = abs(y).^2; % Energy of received signal over N samples
 Test_ed = (1/N(num)).*sum(energy)/(Ps/(snr)); % Test Statistic for the Energy Detection
 thresh_ed = (sqrt(2).*qfuncinv(pfa_ed(m))./sqrt(N(num)))+ 1; % Theoretical value of threshold
 
 pfa_mme = (0.4642 - pfa_ed(m))/(1 - pfa_ed(m));
 part1 = ((sqrt(N(num))+sqrt(M*L))/(sqrt(N(num))-sqrt(M*L)))^2;
 part2 = invtw(1-pfa_mme)*(sqrt(N(num))+sqrt(M*L))^(-2/3);
 part3 = 1 + part2/((N(num)*M*L)^(1/6));
 thresh_mme = part1*part3;  % Theoretical value of threshold
 
 if(Test_ed >= thresh_ed)  % Check whether the received energy is greater than threshold, if so, increment Pd (Probability of detection) detctor by 1
     
     detect = detect+1;
     
 else %Moving to MME
     
    y_mid = [y zeros(1,L)];
    yy = transpose(y_mid);

    summ = zeros(M*L,M*L);
    for i=1:N(num)+1
        
        y_h =  flip(yy(i:L+i-1,:));
        y_hat = reshape(y_h,[M*L,1]);
        summ = summ + y_hat*y_hat';

    end
    
    R = summ/N(num);  %Statistical Covariance Matrix
    lamda = eig(R); %Eigenvalues of R
    Test_mme = max(lamda)/min(lamda);% Test Statistic for MME
    
    
    if(Test_mme >= thresh_mme)
        
        detect = detect+1;
        
    end
     
 end
 
end
Pd(num,m) = detect/k; %Computing the Probability of Detection
end
end

Pd_tot = Pd.^3;

figure;
plot(pfa_ed,Pd_tot(6:11,:),'-o')
grid on
title("Pd Vs Pf")
xlabel("Probability of False Alarm for Energy Detection Stage")
ylabel("Cooperative Probability of Detction Pd")
legend('N=30000','N=32000','N=34000','N=36000','N=38000','N=40000');



surf(pfa_ed,N(6:11),Pd_tot(6:11,:))
xlabel("Probability of False Alarm for Energy Detection Stage")
ylabel("Number of Samples")
%%
clc
clear 
close all

L=32; %Smoothing Factor
M=1; %Number of Antennas
P = M; 
pfa_ed = 0:0.05:0.45; %Probability of False Alarm
N = 6e3:1e3:1*1e4; %Number of Samples
snr_dB = -21; % SNR in decibels
snr = 10^(snr_dB./10); % Linear Value of SNR
%% Simulation to plot Probability of Detection (Pd) vs. Signal to Noise Ratio (SNR) 
%Changing Number of samples

for num=1:length(N)
    
for m = 1:length(pfa_ed)
    
    detect = 0;
for k=1:2000 % Number of Monte Carlo Simulations
    
    n=0:N(num)-1;
    S = 1.*sin(2*pi*1/N(num)*n); %A sinusoidal signal
    Ps = sum(abs(S).^2)/length(S); %Power of the signal
    y = awgn(S,snr_dB, 'measured'); %Adding the effect of AWGN

 energy = abs(y).^2; % Energy of received signal over N samples
 Test_ed = (1/N(num)).*sum(energy)/(Ps/(snr)); % Test Statistic for the Energy Detection
 thresh_ed = (sqrt(2).*qfuncinv(pfa_ed(m))./sqrt(N(num)))+ 1; % Theoretical value of threshold
 
 pfa_mme = (0.4642 - pfa_ed(m))/(1 - pfa_ed(m));
 part1 = ((sqrt(N(num))+sqrt(M*L))/(sqrt(N(num))-sqrt(M*L)))^2;
 part2 = invtw(1-pfa_mme)*(sqrt(N(num))+sqrt(M*L))^(-2/3);
 part3 = 1 + part2/((N(num)*M*L)^(1/6));
 thresh_mme = part1*part3;  % Theoretical value of threshold
 
 if(Test_ed >= thresh_ed)  % Check whether the received energy is greater than threshold, if so, increment Pd (Probability of detection) detctor by 1
     
     detect = detect+1;
     
 else %Moving to MME
     
    y_mid = [y zeros(1,L)];
    yy = transpose(y_mid);

    summ = zeros(M*L,M*L);
    for i=1:N(num)+1
        
        y_h =  flip(yy(i:L+i-1,:));
        y_hat = reshape(y_h,[M*L,1]);
        summ = summ + y_hat*y_hat';

    end
    
    R = summ/N(num);  %Statistical Covariance Matrix
    lamda = eig(R); %Eigenvalues of R
    Test_mme = max(lamda)/min(lamda);% Test Statistic for MME
    
    
    if(Test_mme >= thresh_mme)
        
        detect = detect+1;
        
    end
     
 end
 
end
Pd(num,m) = detect/k; %Computing the Probability of Detection
end
end

Pd_tot = Pd.^3;

figure;
plot(pfa_ed,Pd_tot,'-o')
grid on
title("Pd Vs Pf")
xlabel("Probability of False Alarm for Energy Detection Stage")
ylabel("Cooperative Probability of Detction Pd")
legend('N=6000','N=7000','N=8000','N=9000','N=10000');


figure;
surf(pfa_ed,N,Pd_tot)
xlabel("Probability of False Alarm for Energy Detection Stage")
ylabel("Number of Samples")

%%
%%%%%%% Majority RULE EXPERIMENTS
clc
clear


zfun = @(Pfa1,Pfa2) Pfa1 + Pfa2 - Pfa1.*Pfa2 ;
zhandle = fcontour(zfun);

zhandle.LevelList = 0.1958;
xlabel Pfa_-ed
ylabel Pfa_-mme
title('Pfa_-ed +(1-Pfa_-ed)*Pfa_-mme')
grid on
zhandle.YRange = [0,0.5];
zhandle.XRange = [0,0.5];
axis equal

%%
clc
clear 
close all

L=8; %Smoothing Factor
M=1; %Number of Antennas
P = M; 
pfa_ed = 0:0.01:0.19; %Probability of False Alarm
pfa_mme = 0:0.01:0.19; %Probability of False Alarm
N = 1e4:1000:2*1e4; %Number of Samples
snr_dB = -21; % SNR in decibels
snr = 10^(snr_dB./10); % Linear Value of SNR
%% Simulation to plot Probability of Detection (Pd) vs. Signal to Noise Ratio (SNR) 
%Changing Number of samples

for num=1:length(N)
    
for m = 1:length(pfa_ed)
    
    detect = 0;
for k=1:2000 % Number of Monte Carlo Simulations
    
    n=0:N(num)-1;
    S = 1.*sin(2*pi*1/N(num)*n); %A sinusoidal signal
    Ps = sum(abs(S).^2)/length(S); %Power of the signal
    y = awgn(S,snr_dB, 'measured'); %Adding the effect of AWGN

 energy = abs(y).^2; % Energy of received signal over N samples
 Test_ed = (1/N(num)).*sum(energy)/(Ps/(snr)); % Test Statistic for the Energy Detection
 thresh_ed = (sqrt(2).*qfuncinv(pfa_ed(m))./sqrt(N(num)))+ 1; % Theoretical value of threshold
 
 pfa_mme = (0.4642 - pfa_ed(m))/(1 - pfa_ed(m));
 part1 = ((sqrt(N(num))+sqrt(M*L))/(sqrt(N(num))-sqrt(M*L)))^2;
 part2 = invtw(1-pfa_mme)*(sqrt(N(num))+sqrt(M*L))^(-2/3);
 part3 = 1 + part2/((N(num)*M*L)^(1/6));
 thresh_mme = part1*part3;  % Theoretical value of threshold
 
 if(Test_ed >= thresh_ed)  % Check whether the received energy is greater than threshold, if so, increment Pd (Probability of detection) detctor by 1
     
     detect = detect+1;
     
 else %Moving to MME
     
    y_mid = [y zeros(1,L)];
    yy = transpose(y_mid);

    summ = zeros(M*L,M*L);
    for i=1:N(num)+1
        
        y_h =  flip(yy(i:L+i-1,:));
        y_hat = reshape(y_h,[M*L,1]);
        summ = summ + y_hat*y_hat';

    end
    
    R = summ/N(num);  %Statistical Covariance Matrix
    lamda = eig(R); %Eigenvalues of R
    Test_mme = max(lamda)/min(lamda);% Test Statistic for MME
    
    
    if(Test_mme >= thresh_mme)
        
        detect = detect+1;
        
    end
     
 end
 
end
Pd(num,m) = detect/k; %Computing the Probability of Detection
end
end


Pd_tot = -2.*(Pd.^3) + 3*(Pd.^2);

figure;
plot(pfa_ed,Pd_tot(1:2:end,:),'-o')
grid on
xlim([0,0.19])
title("Pd Vs Pf")
xlabel("Probability of False Alarm for Energy Detection Stage")
ylabel("Cooperative Probability of Detction Pd")
legend('N=10000','N=12000','N=14000','N=16000','N=18000','N=20000')



surf(pfa_ed,N,Pd_tot)
xlabel("Probability of False Alarm for Energy Detection Stage")
ylabel("Number of Samples")

%%
clc
clear 
close all

L=32; %Smoothing Factor
M=1; %Number of Antennas
P = M; 
pfa_ed = 0:0.01:0.19; %Probability of False Alarm
%pfa_mme = 0:0.01:0.19; %Probability of False Alarm
N = 3000:500:5000;%:1000:2*1e4; %Number of Samples
snr_dB = -21; % SNR in decibels
snr = 10^(snr_dB./10); % Linear Value of SNR
%% Simulation to plot Probability of Detection (Pd) vs. Signal to Noise Ratio (SNR) 
%Changing Number of samples

for num=1:length(N)
    
for m = 1:length(pfa_ed)
    
    detect = 0;
for k=1:2000 % Number of Monte Carlo Simulations
    
    n=0:N(num)-1;
    S = 1.*sin(2*pi*1/N(num)*n); %A sinusoidal signal
    Ps = sum(abs(S).^2)/length(S); %Power of the signal
    y = awgn(S,snr_dB, 'measured'); %Adding the effect of AWGN

 energy = abs(y).^2; % Energy of received signal over N samples
 Test_ed = (1/N(num)).*sum(energy)/(Ps/(snr)); % Test Statistic for the Energy Detection
 thresh_ed = (sqrt(2).*qfuncinv(pfa_ed(m))./sqrt(N(num)))+ 1; % Theoretical value of threshold
 
 pfa_mme = (0.4642 - pfa_ed(m))/(1 - pfa_ed(m));
 part1 = ((sqrt(N(num))+sqrt(M*L))/(sqrt(N(num))-sqrt(M*L)))^2;
 part2 = invtw(1-pfa_mme)*(sqrt(N(num))+sqrt(M*L))^(-2/3);
 part3 = 1 + part2/((N(num)*M*L)^(1/6));
 thresh_mme = part1*part3;  % Theoretical value of threshold
 
 if(Test_ed >= thresh_ed)  % Check whether the received energy is greater than threshold, if so, increment Pd (Probability of detection) detctor by 1
     
     detect = detect+1;
     
 else %Moving to MME
     
    y_mid = [y zeros(1,L)];
    yy = transpose(y_mid);

    summ = zeros(M*L,M*L);
    for i=1:N(num)+1
        
        y_h =  flip(yy(i:L+i-1,:));
        y_hat = reshape(y_h,[M*L,1]);
        summ = summ + y_hat*y_hat';

    end
    
    R = summ/N(num);  %Statistical Covariance Matrix
    lamda = eig(R); %Eigenvalues of R
    Test_mme = max(lamda)/min(lamda);% Test Statistic for MME
    
    
    if(Test_mme >= thresh_mme)
        
        detect = detect+1;
        
    end
     
 end
 
end
Pd2(num,m) = detect/k; %Computing the Probability of Detection
end
end


Pd_tot2 = -2.*(Pd2.^3) + 3*(Pd2.^2);


figure;
plot(pfa_ed,Pd_tot,'-o')
grid on
xlim([0,0.19])
title("Pd Vs Pf")
xlabel("Probability of False Alarm for Energy Detection Stage")
ylabel("Cooperative Probability of Detction Pd")
legend('N=3000','N=3500','N=4000','N=4500','N=5000')


figure;
surf(pfa_ed,N,Pd_tot)
xlabel("Probability of False Alarm for Energy Detection Stage")
ylabel("Number of Samples")
