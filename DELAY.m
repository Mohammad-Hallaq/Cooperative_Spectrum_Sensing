clc
close all
clear 
N = 1000;%:1000:1e4; %Number of Samples
snr_dB = -25:5; % SNR in decibels
snr = 10.^(snr_dB./10); % Linear Value of SNR
Pf = 0.1; %Probability of False Alarm
m = 5;
detect = 0;
M = 1;
%% MEASURING DELAY OF ED

for k=1:5000 % Number of Monte Carlo Simulations
    
    n=0:N-1;
   for ind=1:M
        
        S(ind,:) = 1.*sin(2*pi*1/N*n); %A sinusoidal signal
        Ps(ind) = sum(abs(S(ind,:)).^2)/length(S(ind,:)); %Power of the signal
        y(ind,:) = awgn(S(ind,:),snr_dB(m), 'measured'); %Adding the effect of AWGN
       
    end

    tic;
 energy = abs(y).^2; % Energy of received signal over N samples
 Test_ed = (1/(M*N)).*sum(energy)./(Ps(1)/(snr(m))); % Test Statistic for the Energy Detection
 thresh_ed = (sqrt(2).*qfuncinv(0.1)./sqrt(M*N))+ 1; % Theoretical value of threshold
 
 if(Test_ed >= thresh_ed)  % Check whether the received energy is greater than threshold, if so, increment Pd (Probability of detection) detect by 1
     detect = 1;
 else
     detect = 0;
 end
 
 end_time(k) = toc;
end
Pd = detect/k; %Computing the Probability of Detection

mean(end_time)
%%
clc
clear 
close all

L=8; %Smoothing Factor
M=6; %Number of Antennas
P = M; 
pf=0.1; %Probability of False Alarm
N = 1000; %Number of Samples
snr_dB = -25:5; % SNR in decibels
snr = 10.^(snr_dB./10); % Linear Value of SNR
m = 5;
detect = 0;
%% MEASURING DELAY OF MME

for k=1:5000 % Number of Monte Carlo Simulations
    
     n=0:N-1;
    for ind=1:M
        
        S(ind,:) = 1.*sin(2*pi*1/N*n); %A sinusoidal signal
        y(ind,:) = awgn(S(ind,:),snr_dB(m), 'measured'); %Adding the effect of AWGN
    end
    
    tic;
    
    y_mid = horzcat(y ,zeros(M,L));
    yy = transpose(y_mid);
    summ = zeros(M*L,M*L);
     for i=1:N+1
        
        y_h =  flip(yy(i:L+i-1,:));
        y_hat = reshape(y_h,[M*L,1]);
        summ = summ + y_hat*y_hat';

     end

    R = summ/N;  %Statistical Covariance Matrix
    lamda = eig(R); %Eigenvalues of R
    Test_mme = max(lamda)/min(lamda); % Test Statistic for MME
    
 part1 = ((sqrt(N)+sqrt(M*L))/(sqrt(N)-sqrt(M*L)))^2;
 part2 = 0.45*(sqrt(N)+sqrt(M*L))^(-2/3);
 part3 = 1 + part2/((N*M*L)^(1/6));
 thresh = part1*part3;  % Theoretical value of threshold
 
 if(Test_mme >= thresh)  % Check whether the received energy is greater than threshold, if so, increment Pd (Probability of detection) detect by 1
     detect = 1;
 else
     detect = 0;
 end
 
 end_time(k) = toc;
end
Pd = detect/k; %Computing the Probability of Detection

mean(end_time)
%%
clc
clear 
close all

L=32; %Smoothing Factor
M=1; %Number of Antennas
P = M; 
pf=0.1; %Probability of False Alarm
N = 4500; %Number of Samples
snr_dB = -25:5; % SNR in decibels
snr = 10.^(snr_dB./10); % Linear Value of SNR
m = 5;
detect = 0;
num =1 ;
%% MEASURING DELAY OF ED_MME

for k=1:5000 % Number of Monte Carlo Simulations
    
    n=0:N-1;
    for ind=1:M(num)
        
        S(ind,:) = 1.*sin(2*pi*1/N*n); %A sinusoidal signal
        y(ind,:) = awgn(S(ind,:),snr_dB(m), 'measured'); %Adding the effect of AWGN
    end
 
 tic;
 
 Ps = sum(abs(S(1,:).^2))/N; %Power of the signal
 energy = abs(y).^2; % Energy of received signal over N samples
 Test_ed = (1/(M(num)*N)).* sum(sum(energy))/(Ps(1)/(snr(m))); % Test Statistic for the Energy Detection
 thresh_ed = (sqrt(2).*qfuncinv(0.14)./sqrt(M(num)*N))+ 1; % Theoretical value of threshold
 
 
 if(Test_ed >= thresh_ed)  % Check whether the received energy is greater than threshold, if so, increment Pd (Probability of detection) detctor by 1
     
     detect = 1;
     end_time(k) = toc;
     
 else %Moving to MME
   
    y_mid = horzcat(y ,zeros(M(num),L));
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
    part1 = ((sqrt(N)+sqrt(M*L))/(sqrt(N)-sqrt(M*L)))^2;
    part2 = invtw(1-0.064)*(sqrt(N)+sqrt(M*L))^(-2/3);
    part3 = 1 + part2/((N*M*L)^(1/6));
    thresh_mme = part1*part3;  % Theoretical value of threshold
    
    if(Test_mme >= thresh_mme)
        
        detect = 1;
    else
        
        detect = 0;
    end
 
    end_time(k) = toc;
        
  end
     
 end
 
Pd = detect/k; %Computing the Probability of Detection
mean(end_time)
