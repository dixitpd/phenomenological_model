function [o] = JSD(x1,x2)

m  = 0.5*(x1+x2);
t1 = nansum((x1.*log(x1./m))');
t2 = nansum((x2.*log(x2./m))');

o = t1+t2;


end