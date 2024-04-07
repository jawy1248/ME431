clear; clc; close all;
format compact;
p = zeros(1,3);
p(2) = 1;

%% P9-11
if p(1)
    A = [1, 2; 0, 1];
    B = [1; -1];
    C = [0, 3];
    
    Cab = [B, A*B]
    
    syms s;
    
    eq1 = (s+2-3*(1j)) * (s+2+3*(1j));
    poles = simplify(eq1)
    
    charEq = det(s*eye(2) - A)
end
%% P12-16
if p(2)
    % P12
    A = [-1, 4; 1, 0];
    B = [1; 0];
    C = [1, -1];
    
    Cab = [B, A*B]
    
    syms s;
    
    eq1 = (s+1-2*(1j)) * (s+1+2*(1j));
    poles = simplify(eq1)
    
    charEq = det(s*eye(2) - A)
    
    % P13
    K = [2, 6];
    Aa = [1, 1; 0, 1];
    
    cl_Poles = eig(A-B*K)

    % P14
    K = [2, 5]
    Kr = -1/(C*inv(A-B*K)*B)

    % P15
    Obs = [C; C*A]
end

%% P17
if p(3)
    syms s;
    P = -1/(s+5)
end
