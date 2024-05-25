function [R_OSE, varnt1] = OSE_Estimation_new(X,D)

% Description: Optimal Subspace Averaging technique by Vaccaro
% et al.
%
% Input parameters:
% X : NxM matriz containing the snapshots (N antennas, M snapshots)
% D : Presumed dimension of the signal subspace (needed by OSE)
%
% Output parameters:
% R_OSE: NxN estimate of the covariance matrix

[N,M] = size(X);
R_hat = (X*X')/M;


%% Calculation of Z
[U, SS , ~] = svd(X);
U1 = U(:,1:D);
U2 = U(:,D+1:N);
U1up = U1(1:N-1,:);
U1down = U1(2:N,:);
U2up = U2(1:N-1,:);
U2down = U2(2:N,:);
U1plus = (inv(U1down'*U1down))*U1down';
S1 = SS(1:D, 1:D); S1 = S1 + 1e-10*eye(D); % TO AVOID NUMERICAL ISSUE
Pper1 = eye(N-1)-(U1down*U1plus);
W = orth(Pper1);
H = (kron(inv(S1),(W'*U2up)))-kron(((inv(S1))*U1plus*U1up).',(W'*U2down));
r=rank(H);
[U3, S2, V2] = svd(H);
U4 = U3(:,1:r);
S3 = S2(1:r, 1:r);
V3 = V2(:,1:r);
rls = (W'*U1up);
Z = V3*(inv(S3))*U4'*(rls(:));
Zhat = reshape(Z, N-D, D);

%% OSE Subspace Estimation
X1 = U1 - (U2*Zhat*(inv(S1)));
X1 = orth(X1);
P = X1*(inv(X1'*X1))*X1';
Pper = eye(N)-P;

%% Estimation of noise variance
lambda = sort(real(svd(R_hat)),'descend');
lambda = lambda(lambda>1e-6);  % we eliminate very small eigenvalues
                               % (useful when the number of snapshots is smaller than the number of antennas)
lambda = lambda(D+1:end);
varnt1 = ((M)/(M-1))* mean(lambda);


%% R_OSE estimation
R_OSE = (P*R_hat*P) + (varnt1 * Pper) ;  % original OSE

