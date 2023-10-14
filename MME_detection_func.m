function [Pd] = MME_detection_func(Pfa, N, M, L, SNR_dB)


Pd = zeros(1,length(SNR_dB));
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
        y(ind,:) = awgn(S(ind,:),SNR_dB(m), 'measured'); %Adding the effect of AWGN
        y_mid(ind,:) = [y(ind,:) zeros(1,L)];
    end
    
    
    yy = transpose(y_mid);
    
    summ = zeros(M(num)*L,M(num)*L);
    for i=1:N+1
        
        y_h =  flip(yy(i:L+i-1,:));
        y_hat = reshape(y_h,[M(num)*L,1]);
        summ = summ + y_hat*y_hat';
        
    end

    R = summ/N;  %Statistical Covariance Matrix
    lamda = eig(R); %Eigenvalues of R
    Test_mme = max(lamda)/min(lamda); % Test Statistic for MME
    
 part1 = ((sqrt(N)+sqrt(M(num)*L))/(sqrt(N)-sqrt(M(num)*L)))^2;
 part2 = invtw(1-Pfa)*(sqrt(N)+sqrt(M(num)*L))^(-2/3);
 part3 = 1 + part2/((N*M(num)*L)^(1/6));
 thresh = part1*part3;  % Theoretical value of threshold
 
 if(Test_mme >= thresh)  % Check whether the received energy is greater than threshold, if so, increment Pd (Probability of detection) detect by 1
     detect = detect+1;
 end
end

Pd(num,m) = detect/k; %Computing the Probability of Detection
end

end

end