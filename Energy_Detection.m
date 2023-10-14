% Energy Detection

clc
close all
clear 
N = 1000:1000:1e4;%:1000:1e4; %Number of Samples
snr_dB = -25:5; % SNR in decibels
snr = 10.^(snr_dB./10); % Linear Value of SNR
Pf = 0.1; %Probability of False Alarm
%% Simulation to plot Probability of Detection (Pd) vs. Signal to Noise Ratio (SNR) 
%Cahnging Number of Samples


for num=1:length(N)

for m = 1:length(snr)
    detect = 0;
    
for k=1:2000 % Number of Monte Carlo Simulations
    
    n=0:N(num)-1;
    S = 1.*sin(2*pi*1/N(num)*n); %A sinusoidal signal
    Ps = sum(abs(S).^2)/length(S); %Power of the signal
    %y = awgn(S,snr_dB(m),'measured'); %Adding the effect of AWGN
    y = randn(1,N(num));
 energy = abs(y).^2; % Energy of received signal over N samples
 Test_ed = (1/N(num)).*sum(energy)/((snr(m))); % Test Statistic for the Energy Detection
 thresh_ed = (sqrt(2).*qfuncinv(0.1)./sqrt(N(num)))+ 1; % Theoretical value of threshold
 
 if(Test_ed >= thresh_ed)  % Check whether the received energy is greater than threshold, if so, increment Pd (Probability of detection) detect by 1
     detect = detect+1;
 end
 
end
Pd(num,m) = detect/k; %Computing the Probability of Detection
end
end
figure
plot(snr_dB, Pd,'-o')
grid on
title("Pd Vs SNR")
xlabel("Signal to Noise Ratio SNR (dB)")
ylabel("Probability of Detction Pd")
legend('N=1000','N=2000','N=3000','N=4000','N=5000','N=6000','N=7000','N=8000','N=9000','N=10000');

%% Theroretical ecpression of Probability of Detection; refer above reference.
thresh = (sqrt(2).*qfuncinv(Pf)./sqrt(N))+ 1;
Pd_the = qfunc(((thresh - (snr + 1)).*sqrt(N))./(sqrt(2).*sqrt(2*snr + 1)));
plot(snr_dB, Pd_the, 'r')
hold on
plot(snr_dB, Pd, 'b-o')
%%

clc
close all
clear 
N = 1000;%:1000:1e4; %Number of Samples
snr_dB = -25:5; % SNR in decibels
snr = 10.^(snr_dB./10); % Linear Value of SNR
Pf = 0.1; %Probability of False Alarm
M = 1:6;
%% Simulation to plot Probability of Detection (Pd) vs. Signal to Noise Ratio (SNR) 
%Changing Number of antenna

for num=1:length(M)

for m = 1:length(snr)
    detect = 0;
    
for k=1:2000 % Number of Monte Carlo Simulations
    
     n=0:N-1;
    for ind=1:M(num)
        
        S(ind,:) = 1.*sin(2*pi*1/N*n); %A sinusoidal signal
        Ps(ind) = sum(abs(S(ind,:)).^2)/length(S(ind,:)); %Power of the signal
        y(ind,:) = awgn(S(ind,:),snr_dB(m), 'measured'); %Adding the effect of AWGN
       
    end

 energy = abs(y).^2; % Energy of received signal over N samples
 Test_ed = (1/(M(num)*N)).* sum(sum(energy))/(Ps(1)/(snr(m))); % Test Statistic for the Energy Detection
 thresh_ed = (sqrt(2).*qfuncinv(0.1)./sqrt(M(num)*N))+ 1; % Theoretical value of threshold
 
 if(Test_ed >= thresh_ed)  % Check whether the received energy is greater than threshold, if so, increment Pd (Probability of detection) detect by 1
     detect = detect+1;
 end
end
Pd(num,m) = detect/k; %Computing the Probability of Detection
end
end

figure
plot(snr_dB, Pd,'-o')
grid on
title("Pd Vs SNR")
xlabel("Signal to Noise Ratio SNR (dB)")
ylabel("Probability of Detction Pd")
legend('M=1','M=2','M=3','M=4','M=5','M=6');
%%
% Energy Detection
clc
close all
clear 
N = 1000;%:1000:1e4; %Number of Samples
snr_dB = -25:5; % SNR in decibels
snr = 10.^(snr_dB./10); % Linear Value of SNR
Pf = 0.1; %Probability of False Alarm
alpha =[0 1 2];

%% Simulation to plot Probability of Detection (Pd) vs. Signal to Noise Ratio (SNR) 
%Cahnging Noise Uncertainty Factor

for num=1:length(alpha)
    
for m = 1:length(snr)
    detect = 0;
    
for k=1:2000 % Number of Monte Carlo Simulations
    
    n=0:N-1;
    S = 1.*sin(2*pi*1/N*n); %A sinusoidal signal
    Ps = sum(abs(S).^2)/length(S); %Power of the signal
    y = awgn(S,snr_dB(m), 'measured'); %Adding the effect of AWGN
    
 sigma(num) = (Ps*10^(alpha(num)/10)/(snr(m)));
 
 energy = abs(y).^2; % Energy of received signal over N samples
 Test_ed = (1/N).*sum(energy)/(sigma(num)); % Test Statistic for the Energy Detection
 thresh_ed = (sqrt(2).*qfuncinv(Pf)./sqrt(N))+ 1; % Theoretical value of threshold
 
 if(Test_ed >= thresh_ed)  % Check whether the received energy is greater than threshold, if so, increment Pd (Probability of detection) detect by 1
     detect = detect+1;
 end
end
Pd(num,m) = detect/k; %Computing the Probability of Detection
end
end
figure
plot(snr_dB, Pd,'-o')
grid on
title("Pd Vs SNR")
xlabel("Signal to Noise Ratio SNR (dB)")
ylabel("Probability of Detction Pd")
legend('ED-0dB','ED-1dB','ED-2dB');

%%
clc
close all
clear 
N = 2000:2000:3e4;%:1000:1e4; %Number of Samples
snr_dB = -21; % SNR in decibels
snr = 10^(snr_dB/10); % Linear Value of SNR
Pf = 0:0.01:1; %Probability of False Alarm
%% Simulation to plot Probability of Detection (Pd) vs. Signal to Noise Ratio (SNR) 
%Cahnging Number of Samples

for num=1:length(N)

for m = 1:length(Pf)
    detect = 0;
    
for k=1:2000 % Number of Monte Carlo Simulations
    
    n=0:N(num)-1;
    S = 1.*sin(2*pi*1/N(num)*n); %A sinusoidal signal
    Ps = sum(abs(S).^2)/length(S); %Power of the signal
    y = awgn(S,snr_dB, 'measured'); %Adding the effect of AWGN
    
 energy = abs(y).^2; % Energy of received signal over N samples
 Test_ed = (1/N(num)).*sum(energy)/(Ps/(snr)); % Test Statistic for the Energy Detection
 thresh_ed = (sqrt(2).*qfuncinv(Pf(m))./sqrt(N(num)))+ 1; % Theoretical value of threshold
 
 if(Test_ed >= thresh_ed)  % Check whether the received energy is greater than threshold, if so, increment Pd (Probability of detection) detect by 1
     detect = detect+1;
 end
 
end
Pd(num,m) = detect/k; %Computing the Probability of Detection
end
end
figure
plot(Pf, Pd,'-o')
grid on
title("Pd Vs SNR")
xlabel("Signal to Noise Ratio SNR (dB)")
ylabel("Probability of Detction Pd")
legend('N=1000','N=2000','N=3000','N=4000','N=5000','N=6000','N=7000','N=8000','N=9000','N=10000');

