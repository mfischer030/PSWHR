
% We can have only two equations in terms of two unknowns: n_r and Rs_r

function [res,Rs_r,n_r]= TwoParameters(x,Ns,k,Tref,q, beta, T, alpha,c1,c2,c3,c4,c5, KT)

Rs_r = x(1);
n_r  = x(2);

res(1)= c3*(2*c4-c5)/(c3-c4*Rs_r)+c4*(c2-2*c3)/(c3-Rs_r*c4)*exp((c5*Rs_r-c2)/n_r) + ((c5*c3-c4*c2)/(c3-Rs_r*c4)-(c2*c4+c5*c3-c5*c2)/n_r)*exp((c3+Rs_r*c4-c2)/n_r);

res(2)= (alpha*(T-Tref)-(beta*(T-Tref)*c4)/(Ns*(c3-c4*Rs_r))+(c1*c4*(2*c3-c2))/(c3-Rs_r*c4)) + (alpha*(T-Tref)*((c3+Rs_r*c4-c2)/n_r-1))*exp((c3+Rs_r*c4-c2)/n_r)...
        + (c4*beta*(T-Tref)/Ns)*((n_r-Rs_r*c4+c3)/(n_r*(c3-Rs_r*c4)))*exp((c3+Rs_r*c4-c2)/n_r) + c1*c4*KT*(c2-2*c3)/(c3-Rs_r*c4)*exp(q*beta*(T-Tref)/(Ns*k*T*n_r)+c2/n_r*(Tref/T-1))...
        + c1*c4*(KT-1)*(2*c3-c2)/(c3-Rs_r*c4)*exp(-c2/n_r);
      
end
