function y = inv_K_N (x,N,K)

inp = 0:0.0001:1;

min = 1;



for i=1:length(inp)
    
    pfa_K_N = 0;
    for j=K:N
    
    pfa_K_N = pfa_K_N + nchoosek(N,j)*(inp(i)^j)*((1-inp(i))^(N-j));
    
    end
    
    if abs(pfa_K_N - x)< min
        
        min = abs(pfa_K_N - x);
        y = inp(i);
        
    
     
    end
   

end



end
