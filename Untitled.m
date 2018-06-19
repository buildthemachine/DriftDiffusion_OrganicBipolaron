clear all
close all

f = 1;
y = vapsolve( y*(1*q*y)^2* (exp(g1_pre*(2*coeff_1*y/Nt)^delta) * exp(g2_pre*(sqrt(1+0.8*(a/esig*coeff_2*fmin)^2)-1)) * y)^2 ...,
        + y*(1*q*y)^2 * exp(g1_pre*(2*coeff_1*y/Nt)^delta) * exp(g2_pre*(sqrt(1+0.8*(a/esig*coeff_2*f)^2)-1)) * y ...,
        * dg2_pre/sqrt(1+0.8*(a/esig*coeff_2*f)^2)*coeff_2*f ...,
        + (1*q*y)*( exp(g1_pre*(2*coeff_1*y/Nt)^delta)*exp(g2_pre*(sqrt(1+0.8*(a/esig*coeff_2*f)^2)-1))*y*f -1) * (1+y*dg1_pre*(2*coeff_1*y/Nt)^(delta-1) ) ...,
        - q*y* (exp(g1_pre*(2*coeff_1*y/Nt)^delta)*exp(g2_pre*(sqrt(1+0.8*(a/esig*coeff_2*f)^2)-1))*y*f -1)^2 == 0, f);    
