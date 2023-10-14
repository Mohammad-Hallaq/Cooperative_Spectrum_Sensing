clc
clear 
close all

L=8; %Smoothing Factor
M=1; %Number of Antennas
P = M; 
pf=0.1; %Probability of False Alarm
N = 1000:1000:1e4; %Number of Samples
snr_dB = -25:5; % SNR in decibels
snr = 10.^(snr_dB./10); % Linear Value of SNR
%% Simulation to plot Probability of Detection (Pd) vs. Signal to Noise Ratio (SNR) 
%Changing Number of samples

for num=1:length(N)
    
for m = 1:length(snr)
    detect = 0;
for k=1:2000 % Number of Monte Carlo Simulations
    
    n=0:N(num)-1;
    S = 1.*sin(2*pi*1/N(num)*n); %A sinusoidal signal
    Ps = sum(abs(S).^2)/length(S); %Power of the signal
    y = awgn(S,snr_dB(m), 'measured'); %Adding the effect of AWGN

 energy = abs(y).^2; % Energy of received signal over N samples
 Test_ed = (1/N(num)).*sum(energy)/(Ps/(snr(m))); % Test Statistic for the Energy Detection
 thresh_ed = (sqrt(2).*qfuncinv(0.09)./sqrt(N(num)))+ 1; % Theoretical value of threshold
 
 
 part1 = ((sqrt(N(num))+sqrt(M*L))/(sqrt(N(num))-sqrt(M*L)))^2;
 part2 = 2.02*(sqrt(N(num))+sqrt(M*L))^(-2/3);
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

figure;
plot(snr_dB,Pd,'-o')
grid on
title("Pd Vs SNR")
xlabel("Signal to Noise Ratio SNR (dB)")
ylabel("Probability of Detction Pd")
legend('N=1000','N=2000','N=3000','N=4000','N=5000','N=6000','N=7000','N=8000','N=9000','N=10000');
%%
clc
clear 
close all

L=8; %Smoothing Factor
M = 1:6; %Number of Antennas
P = M; 
pf=0.1; %Probability of False Alarm
N = 1000; %Number of Samples
snr_dB = -25:5; % SNR in decibels
snr = 10.^(snr_dB./10); % Linear Value of SNR
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
        y_mid(ind,:) = [y(ind,:) zeros(1,L)];
    end

 energy = abs(y).^2; % Energy of received signal over N samples
 Test_ed = (1/(M(num)*N)).* sum(sum(energy))/(Ps(1)/(snr(m))); % Test Statistic for the Energy Detection
 thresh_ed = (sqrt(2).*qfuncinv(0.09)./sqrt(M(num)*N))+ 1; % Theoretical value of threshold
 
 part1 = ((sqrt(N)+sqrt(M(num)*L))/(sqrt(N)-sqrt(M(num)*L)))^2;
 part2 = 2.02*(sqrt(N)+sqrt(M(num)*L))^(-2/3);
 part3 = 1 + part2/((N*M(num)*L)^(1/6));
 thresh_mme = part1*part3;  % Theoretical value of threshold
 
 if(Test_ed >= thresh_ed)  % Check whether the received energy is greater than threshold, if so, increment Pd (Probability of detection) detect by 1
     
     detect = detect+1;
     
 else %Moving to MME
     
    
    yy = transpose(y_mid);

    summ = zeros(M(num)*L,M(num)*L); 
    for i=1:N+1
        
        y_h =  flip(yy(i:L+i-1,:));
        y_hat = reshape(y_h,[M(num)*L,1]);
        summ = summ + y_hat*y_hat';

    end
    
    R = summ/N;  %Statistical Covariance Matrix
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

figure;
plot(snr_dB,Pd,'-o')
grid on
title("Pd Vs SNR")
xlabel("Signal to Noise Ratio SNR (dB)")
ylabel("Probability of Detction Pd")
legend('M=1','M=2','M=3','M=4','M=5','M=6');

%%
clc
clear 
close all

L=2:2:20; %Smoothing Factor
M=1; %Number of Antennas
P = M; 
pf=0.1; %Probability of False Alarm
N = 1000; %Number of Samples
snr_dB = -25:5; % SNR in decibels
snr = 10.^(snr_dB./10); % Linear Value of SNR

%% Simulation to plot Probability of Detection (Pd) vs. Signal to Noise Ratio (SNR) 
%Changing Smoothing Factor

for num=1:length(L)
    
for m = 1:length(snr)
    detect = 0;
for k=1:2000 % Number of Monte Carlo Simulations
    
    n=0:N-1;
    S = 1.*sin(2*pi*1/N*n); %A sinusoidal signal
    Ps = sum(abs(S).^2)/length(S); %Power of the signal
    y = awgn(S,snr_dB(m), 'measured'); %Adding the effect of AWGN

 energy = abs(y).^2; % Energy of received signal over N samples
 Test_ed = (1/N).*sum(energy)/(Ps/(snr(m))); % Test Statistic for the Energy Detection
 thresh_ed = (sqrt(2).*qfuncinv(0.09)./sqrt(N))+ 1; % Theoretical value of threshold
 
 
 part1 = ((sqrt(N)+sqrt(M*L(num)))/(sqrt(N)-sqrt(M*L(num))))^2;
 part2 = 2.02*(sqrt(N)+sqrt(M*L(num)))^(-2/3);
 part3 = 1 + part2/((N*M*L(num))^(1/6));
 thresh_mme = part1*part3;  % Theoretical value of threshold
 
 if(Test_ed >= thresh_ed)  % Check whether the received energy is greater than threshold, if so, increment Pd (Probability of detection) detctor by 1
     
     detect = detect+1;
     
 else %Moving to MME
     
    y_mid = [y zeros(1,L(num))];
    yy = transpose(y_mid);

    summ = zeros(M*L(num),M*L(num));
    for i=1:N+1
      
        y_h =  yy(i:L(num)+i-1,:);
        y_medium = y_h;
         
        for j=1:L(num)
            y_h(j,:)=y_medium(L(num)-j+1,:);
        end
        
        y_hat = reshape(y_h,[M*L(num),1]);
        summ = summ + y_hat*y_hat';

    end

    
    R = summ/N;  %Statistical Covariance Matrix
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

figure;
plot(snr_dB,Pd,'-o')
grid on
title("Pd Vs SNR")
xlabel("Signal to Noise Ratio SNR (dB)")
ylabel("Probability of Detction Pd")
legend('L=2','L=4','L=6','L=8','L=10','L=12','L=14','L=16','L=18','L=20');
%%
clc
clear 
close all

L=8; %Smoothing Factor
M=1; %Number of Antennas
P = M; 
pf=0.1; %Probability of False Alarm
N = 1000; %Number of Samples
snr_dB = -25:5; % SNR in decibels
snr = 10.^(snr_dB./10); % Linear Value of SNR
num =1;
%% Simulation to plot Probability of Detection (Pd) vs. Signal to Noise Ratio (SNR) 
%Changing Number of antennas

for m = 1:length(snr)
    detect = 0;
for k=1:2000 % Number of Monte Carlo Simulations
    
    n=0:N(num)-1;
    S = 1.*sin(2*pi*1/N(num)*n); %A sinusoidal signal
    Ps = sum(abs(S).^2)/length(S); %Power of the signal
    y = awgn(S,snr_dB(m), 'measured'); %Adding the effect of AWGN

 energy = abs(y).^2; % Energy of received signal over N samples
 Test_ed = (1/N(num)).*sum(energy)/(Ps/(snr(m))); % Test Statistic for the Energy Detection
 thresh_ed = (sqrt(2).*qfuncinv(0.09)./sqrt(N(num)))+ 1; % Theoretical value of threshold
 
 
 part1 = ((sqrt(N(num))+sqrt(M*L))/(sqrt(N(num))-sqrt(M*L)))^2;
 part2 = 2.02*(sqrt(N(num))+sqrt(M*L))^(-2/3);
 part3 = 1 + part2/((N(num)*M*L)^(1/6));
 thresh_mme = part1*part3;  % Theoretical value of threshold
 
 if(Test_ed >= thresh_ed)  % Check whether the received energy is greater than threshold, if so, increment Pd (Probability of detection) detctor by 1
     
     detect = detect+1;
     
 else %Moving to MME
     
    y_mid = [y zeros(1,L)];
    yy = transpose(y_mid);

    summ = zeros(M*L,M*L);
    for i=1:N(num)+1
      
        y_h =  yy(i:L+i-1,:);
        y_medium = y_h;
         
        for j=1:L
            y_h(j,:)=y_medium(L-j+1,:);
        end
        
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

%%
clc
clear 
close all

L=8; %Smoothing Factor
M=1; %Number of Antennas
P = M; 
pf=0.1; %Probability of False Alarm
N = 1000; %Number of Samples
snr_dB = -25:5; % SNR in decibels
snr = 10.^(snr_dB./10); % Linear Value of SNR
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
 thresh_ed = (sqrt(2).*qfuncinv(0.09)./sqrt(N))+ 1; % Theoretical value of threshold
 
 
 part1 = ((sqrt(N)+sqrt(M*L))/(sqrt(N)-sqrt(M*L)))^2;
 part2 = 2.02*(sqrt(N)+sqrt(M*L))^(-2/3);
 part3 = 1 + part2/((N*M*L)^(1/6));
 thresh_mme = part1*part3;  % Theoretical value of threshold
 
 if(Test_ed >= thresh_ed)  % Check whether the received energy is greater than threshold, if so, increment Pd (Probability of detection) detctor by 1
     
     detect = detect+1;
     
 else %Moving to MME
     
    y_mid = [y zeros(1,L)];
    yy = transpose(y_mid);

    summ = zeros(M*L,M*L);
    for i=1:N+1
        
        y_h =  flip(yy(i:L+i-1,:));
        y_hat = reshape(y_h,[M*L,1]);
        summ = summ + y_hat*y_hat';

    end
    
    R = summ/N;  %Statistical Covariance Matrix
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

figure;
plot(snr_dB,Pd,'-o')
grid on
title("Pd Vs SNR")
xlabel("Signal to Noise Ratio SNR (dB)")
ylabel("Probability of Detction Pd")
legend('ED-0dB','ED-1dB','ED-2dB');
%%
