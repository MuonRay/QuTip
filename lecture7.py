# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 18:38:40 2022

@author: cosmi
"""


clear;

H = randH(10);

value_primal = inf_norm(H)
value_dual = inf_norm_dual(H)

max(svd(H))


#prime SDP for 
function value = inf_norm(H)
    
    d = max(size(H));
    
    cvx_begin sdp quiet
    
    variable X1(d,d) complex semidefinite
    variable X2(d,d) complex semidefinite
    
    
    maximize trace(H*(X1-X2))
    
    subject to  
          
          trace(X1 + X2) <= 1;
          
    cvx_end
    
    value = cvx_optval;
    
end

## dual SDP for 
function value = inf_norm_dual(H)
    
    d = max(size(H));
    
    cvx_begin sdp quiet
    
    variable t nonnegative
    
    minimize t
    subject to    
    
          H <= t * eye(d);
          
          -t * eye(d) <= H;
          
    cvx_end
    
    value = cvx_optval;
    
end

function H = randH(n)

    H = 2*(randn(n) + i*randn(n)) - (1+i);
    H = H + H';
    
end