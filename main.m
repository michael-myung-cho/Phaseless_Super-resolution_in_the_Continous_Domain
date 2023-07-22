%% Phaseless Super-resolution - 1D case
% 
% Phaseless Atoimc Norm Minimization (PANM) 
% Standard Atomic Norm Minimization (SANM)
% 
% Reference:
% [1] M. Cho, et al, ?Phaseless super-resolution in the continuous domain,? 
% in Proceedings of International Conference on Acoustics, Speech, and Signal Processing (ICASSP), 2017, pp. 3814-3818. 
% [2] G. Tang, et al, ?Compressed sensing off the grid,? 
% IEEE Transactions on Information The- ory, vol. 59, no. 11, pp. 7465?7490, 2013.
% 
% Feb. 18, 2018
% Myung (Michael) Cho 
% michael.myung.cho@gmail.com
%---------------------------------------------------

%% clear 
clear; 
clc;

%% Dimensions of the problem
n = 32;   % signal dimension
ks = 3;   % number of frequencies
ms = 30;   % number of measurements


%% initialization
res =[];
nTrial=1;

for k = ks
    for m = ms
        m
        nSuccPANM = 0; 
        nSuccSANM = 0; 
        for iTrial=1:nTrial
            iTrial
            %% Create a signal with s different frequencies 
            [x0, f0, V0, c0] = sigGen(n,k);
            f0
            
            %% Select samples uniformly at random in the interval
            unif_rand_samples = randperm(n);
            T = (1:n)';         % total time sample
            Omega = T(unif_rand_samples(1:m));
            Omega = sort(Omega);        % randomly chosen time sample
            Omega=(1:m)'; 
            
            %% Phaseless Atomic Norm Minimization         
            [ED_F_PANM]  = PANM(x0, Omega, n, k, f0); 
            if (ED_F_PANM < 10^-3) 
                 nSuccPANM = nSuccPANM + 1;
            end 
            
            %% Standard Atomic Norm Minimization
            [ED_F_SANM]  = SANM_cvx(x0, Omega, n, k, f0); 
            if (ED_F_SANM < 10^-3) 
                 nSuccSANM = nSuccSANM + 1;
            end                       
             
        end
        resBuf = [n,k,m, nSuccPANM, nSuccSANM ];
    end
end










