function [y] = invtw(x)

inp = -4:0.05:2.03;
ftw_x = tracywidom_appx(inp,1);
min= 1;

for i=1:length(inp)
    
    
    if abs(tracywidom_appx(inp(i),1)-x)< min
        
        min = abs(tracywidom_appx(inp(i),1)-x);
        y = inp(i);
        
    end
   

end




end