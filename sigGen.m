function [x0, f0, V0, c0] = sigGen(n,k)
%  Outputs:
%  X0: Original Signal
%  f0: Original frequencies 
%  c0: Original frequency coefficient
%  V0: Frequncy (or Time) atom matrix

% Pick s freqs with the specified spacing type in the interval [0,1]
f0=rand(k,1); 

% Initialize the frequency signal
x0=zeros(n,1);
V0=exp(1i*2*pi*kron((0:n-1)',f0'));
c0=rand(k,1)+1i*rand(k,1);
x0=V0*c0;

end