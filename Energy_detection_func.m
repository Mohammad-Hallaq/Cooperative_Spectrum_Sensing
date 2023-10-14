function [Pd] = Energy_detection_func(Pfa, N, M, SNR_dB)

SNR = 10.^(SNR_dB./10);
for num=1:length(M)

for m = 1:length(SNR)
    detect = 0;
    
for k=1:2000 % Number of Monte Carlo Simulations
    
     n=0:N-1;
     S = zeros(M(num),N);
     y = zeros(M(num),N);
     
    for ind=1:M(num)
        
        S(ind,:) = 1.*sin(2*pi*1/N*n); %A sinusoidal signal
        Ps = sum(abs(S(ind,:)).^2)/length(S(ind,:)); %Power of the signal
        y(ind,:) = awgn(S(ind,:),SNR_dB(m), 'measured'); %Adding the effect of AWGN
       
    end

 energy = abs(y).^2; % Energy of received signal over N samples
 Test_ed = (1/(M(num)*N)).* sum(sum(energy))/(Ps(1)/(SNR(m))); % Test Statistic for the Energy Detection
 thresh_ed = (sqrt(2).*qfuncinv(Pfa)./sqrt(M(num)*N))+ 1; % Theoretical value of threshold
 
 if(Test_ed >= thresh_ed)  % Check whether the received energy is greater than threshold, if so, increment Pd (Probability of detection) detect by 1
     detect = detect+1;
 end
end
Pd(num,m) = detect/k; %Computing the Probability of Detection
end
end

end