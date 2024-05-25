function [Khat, Rsimc] = Order_Estimation_SIMC(Xd,Kmax,mu,itrmax)

%code for paper
% V. Garg, A. Pagès-Zamora, and I. Santamaria,“Order estimation with missing data for massive MIMO systems”, Submitted to the IEEE Signal Processing Letters, 2021.

%Inputs:
%% Xd = Matrix with missing entries
%% Kmax = an overestimation of the rank

%outputs
%% Khat = estimated order/rank
%% Rsimc = Matrix with all the entries after applying SIMC with rank Khat.


[M, ~] = size(Xd);

for p = 1:Kmax
    %% alpha
         Xd1 = Xd(1:M-1,:); Xd2 = Xd(2:M,:);
        [Ud1, ~] = svds(Xd1,p); [Ud2, ~] = svds(Xd2,p);
        alpha = norm(Ud1*Ud1'-Ud2*Ud2','fro')*M/20;
        %% SIMC (with OSE as explained in the paper)
        [Rsimc] = Shift_Invariance_Matrix_Completion(Xd,p,mu,alpha,itrmax);
        %% subspaces
        [Uup, ~] = svds(Rsimc(1:M-1,1:M-1),p);
        [Udown, ~] = svds(Rsimc(2:M,2:M),p);
        %% subspace Distance
        Dsc1(p) = norm(Uup*Uup'-Udown*Udown','fro')/p;

end
[~, Khat] = min(Dsc1);
