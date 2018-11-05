% M = E(:,:,:);
% M = permute(M,[3 1 2]);

M = reshape(permute(E,[3 2 1]), [numel(t) numel(T)]);
M(:,T(:)==0 & vec(eye(size(T))==0) ) = []; % remove out of diagonal kinetics if not active
% M(:,T(:)==0 ) = []; % remove out of diagonal kinetics if not active

plot(M)
legend({'1','2','3','4','5','6'})

%%

Q = M(:,2);
N = M(:,[1 3 6]);

Np=pinv(N);

figure(5)
plot(N*(Np*Q)-Q)
% hold on
% plot(Q)
% hold off