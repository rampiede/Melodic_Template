function [epsilon,F_aprox] = Binary_search_epsilon(F,p )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    [epsilon_cand,F_aprox,~] = Epsilon_candidate(F, p );
    eps_cand = epsilon_cand;
    n = length(epsilon_cand);
    tol = 10^-10;
    while n ~= 1
        i = ceil(n/2);
        epsilon = epsilon_cand(i);
        hold off
        if Decision_Problem(F_aprox, epsilon, p)
            epsilon_cand(i+1:n) = [];
        else
            epsilon_cand(1:i) = [];
        end
        n = length(epsilon_cand);       
    end  
    epsilon = epsilon_cand;
    t = find(eps_cand == epsilon);
    epsilon2 = eps_cand(t-1);%eps_cand(1) == 0, y nunca va a ser epsilon = 0
    len = epsilon - epsilon2;
    while len > tol
        len = len/2;
        epsilon3 = epsilon2 + len;
        if Decision_Problem(F_aprox,epsilon3,p)
            epsilon = epsilon3;
        else
            epsilon2 = epsilon3;
        end
    end
            
end

