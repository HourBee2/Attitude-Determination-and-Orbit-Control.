wn = 1;
damp = 1;

% Underdamped system
sys = tf([0 0 1],[0 1 0]) * tf([0 0 wn^2],[1 2*damp*wn wn^2]); 
rlocus(sys)

clc; clear all
%% Problem 1
I = 1;    % kg m²
Mp = 0.2;
ts = 60;  % s

syms theta u
eq = diff(diff(theta)) * I - u ;

% Parámetros deseados
syms dzeta 
eq = Mp == exp(-pi*dzeta/(1-dzeta)^0.5);
[sol] = vpasolve(eq,dzeta);

%dzeta = sol.dzeta

%Definir tf genérica
s = syms2tf(eq);

%Bucle cerrado 

s = s/(1+s);

%Encontrar polos y zeros
z = zero(s);
p = pole(s);

pzplot(s)
grid

%Definir PD
Kp = 1;
Ki = 1;
Kd = 1;
C = Kp + Kd*s

function[TF] = syms2tf(eq)
% Get transfer function from equation
% 
%   input eq (sym): symbolic equation
% 
%   output TF (tf): transfer function
%

  G = laplace(eq);             %Get Laplace transform of Eq
  [symNum,symDen] = numden(G); %Get num and den of Symbolic TF
  TFnum = sym2poly(symNum);    %Convert Symbolic num to polynomial
  TFden = sym2poly(symDen);    %Convert Symbolic den to polynomial
  TF =tf(TFnum,TFden);
  
end