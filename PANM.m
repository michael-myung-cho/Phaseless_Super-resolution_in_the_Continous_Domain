function [EDFreq] = PANM(x0, M, n, k, f0)
%% Phaseless Atoimc Norm Minimization (PANM)
% 
% x0: grouth truth
% M: number of measurement
% n: signal dimension
% k: number of frequencies
% f0: ground truth freqeuncey
%
% by Myung (Michael) Cho
%--------------------------------------------

Xabs=(abs(x0)).^2;
Mc = (1:n)';
Mc(M) = [];
m=max(size(M));

Ms=[M(2:end);M(1)];
XplusXsAbs=(abs(x0(M)+x0(Ms))).^2;
XminusiXsAbs=(abs(x0(M)-1i*x0(Ms))).^2;


% solve the dual problem
cvx_quiet true
cvx_precision default
%cvx_solver sdpt3

cvx_begin sdp
  variable Q(n,n) hermitian;
  variable u(n) complex;
  dual variable M1;
  dual variable M2;
  
  Tu = toeplitz(u);
  
  minimize real(trace(Tu))
  subject to 
        Tu-Q>=0: M1;
        Q>=0: M2;
        for ii=1:m
            Q(M(ii),M(ii))==Xabs(M(ii));
            Q(M(ii),M(ii))+Q(Ms(ii),Ms(ii))+Q(M(ii),Ms(ii))+Q(Ms(ii),M(ii))==XplusXsAbs(ii);
            Q(M(ii),M(ii))+Q(Ms(ii),Ms(ii))-1i*Q(Ms(ii),M(ii))+1i*Q(M(ii),Ms(ii))==XminusiXsAbs(ii);
        end   
cvx_end  

[fEst, amp] = VanDec(u);

amp = amp * 2;
fEst=sort(fEst,'descend');

[U D V]=svds(Q,1);
xEst=U*sqrt(D);

%% Global Phase Calibration 
cvx_begin 
    variable e complex
    minimize (((e.*xEst - x0)'*(e.*xEst - x0))/(x0'*x0))
cvx_end
norm(e,2)

xHat=e*xEst;

%% estimate the frequencies using esprit
res_err = norm(xHat-x0)/norm(x0);
fprintf('PANM: Relative Error  = %f \n', res_err);
fprintf('---------------------------------------------------------\n');

%% display Mean Square Error (MSE) in frequency
fprintf('\n\n');
nfEst = max(size(fEst));
fErr = 0;
if nfEst >= k
    freqC = [f0;min(f0)+1;max(f0)-1];
    for i=1:nfEst
    % freq error
    [res, ind] = sort(abs(freqC - fEst(i)),'ascend');
    fErr = fErr + res(1)^2;
    end
else 
    fEstC = [fEst;min(f0)+1;max(f0)-1]
    for i=1:k    
        [res, ~] = sort(abs(f0(i) - fEstC),'ascend');
        fErr = fErr + res(1)^2;
    end
end
EDFreq = sqrt(fErr);
fprintf('PANM: Euclidean Distance of freq = %f\n', EDFreq);
fprintf('---------------------------------------------------------\n');

end