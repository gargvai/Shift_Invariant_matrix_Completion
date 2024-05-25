function Rsimc  = Shift_Invariance_Matrix_Completion(Xd,p,mu,alpha,itrmax)
%% Shift Invariance Matrix Completion
%% for the papers
%% V. Garg, P. Giménez-Febrer, A. Pagès-Zamora, and I. Santamaria,“DOA estimation via shift-invariant matrix completion”. Signal Processing, volume 183, 2021.
%% V. Garg, A. Pagès-Zamora, and I. Santamaria,“Order estimation with missing data for massive MIMO systems”, Submitted to the IEEE Signal Processing Letters, 2021.

%% Xd = Matrix with missing entries
%% Rsimc = Output matrix with all the entries after applying SIMC method
%% p = low rank of the matrix

[M N] = size(Xd);
R0 = Xd; # Data matrix with missing entries
Rm = real(R0./R0); Rm = Rm==1;
%% Initializing factors
[UU LL VV] = svds(Xd,p);
W = UU*LL^(1/2); H = VV*LL^(1/2); Q = diag(eye(p))';
%% Initializing Other Parameters
Ratio = 1;
C = 1; %condition
itr = 1;
%% Start Iteration
while C == 1
    for ii = 1:M
        idx = find(real(Rm(ii,:)) == 1);

        if ii == 1
            W(ii,:) = (R0(ii,:)*H + alpha*W(ii+1,:)*diag(Q)) * inv(H(idx,:)'*H(idx,:) + (mu+alpha)*eye(p));
        elseif ii<M
             W(ii,:) = (R0(ii,:)*H + alpha*W(ii-1,:)*diag(Q)'+alpha*W(ii+1,:)*diag(Q)) * inv(H(idx,:)'*H(idx,:) + (mu+alpha)*eye(p) + alpha*diag(Q)*diag(Q)');
        else
            W(ii,:) = (R0(ii,:)*H + alpha*W(ii-1,:)*diag(Q)') * inv(H(idx,:)'*H(idx,:) + mu*eye(p) + alpha*diag(Q)*diag(Q)');
        end
    end
    for jj=1:N
         idx = find(real(Rm(:,jj)) == 1);
         H(jj,:) =(R0(:,jj)'*W)*inv(W(idx,:)'*W(idx,:) + mu*eye(p));
    end
 %% Q
    for pp = 1:p
%%
   AA = 0; BB = 0;
   for ii = 2:M
        AA = AA + W(ii,p)*W(ii,p)';
        BB = BB + W(ii-1,p)*W(ii,p)';
   end
    Q(pp) = (BB) * inv(AA);
    end
%%
  itr = itr+1;
  Ratio_old = Ratio;
  Xmc = W*H';
%% convergence
 R00 = Xmc;
 Ratio = norm((R00.*Rm-R0),'fro')/norm(R0,'fro');
  if abs(Ratio_old-Ratio)/Ratio_old < 1e-10 || itr > itrmax
      C = 0;
  end
  end


Rsimc = OSE_Estimation_new(Xmc,p);

