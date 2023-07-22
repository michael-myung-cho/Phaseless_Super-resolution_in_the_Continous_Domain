function [EDFreq]=SANM_cvx(x0, M, n, k, f0)
%% Standard Atomic Norm Minimization

% solve the dual problem
cvx_quiet true
cvx_precision default

cvx_begin sdp
  variable xHat(n) complex;
  variable u(n) complex;
  variable s complex;
  Q = toeplitz(u);
  
  minimize real(trace(Q)/n + trace(s))
  subject to 
        [s xHat'; xHat Q] >= 0,
        xHat(M)==x0(M);
cvx_end

[fEst, amp] = VanDec(u);

%% estimate the frequencies using esprit
res_err = norm(xHat-x0)/norm(x0);
fprintf('SANM: Relative Error  = %f \n', res_err);
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
fprintf('SANM: Euclidean Distance of freq = %f\n', EDFreq);
fprintf('---------------------------------------------------------\n');

end