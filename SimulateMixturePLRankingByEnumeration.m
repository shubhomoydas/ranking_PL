
%
% Simulate Ranking by Plackett-Luce Model
%
% Assume only 5 items being ranked.
% 1. First, fix v_i's.
% 2. Generate all permutations of the five
%    items (representing each ranking order.)
% 3. Compute the probability for each 
%    permutation and represent as a 120-vector.
% 4. Sample the permutation indexes on the
%    basis of this multinomial distribution.

%%
rng('default'); % MATLAB
%V = unifrnd(0,1,5,1);
V = [
    0.10,0.56,0.30,0.89,0.02;
    0.50,0.10,0.15,0.05,0.20;
    0.25,0.30,0.10,0.10,0.25];
V = diag(1./sum(V,2))*V;
I = 5:-1:1;
PermI = perms(I); % P(1:10,:)
np = size(PermI,1);
nv = size(V,1);
P = zeros(np,nv);
for v=1:nv
    for i=1:np
        tpi = PermI(i,:);
        %tpi = [2,5,3,4,1];
        P(i,v) = computePLProb(tpi,V(v,:));
    end;
end;
%%

rng('default'); % MATLAB
n = 100; % No. samples per mixture component
ni = size(PermI,2); % No. items
nv = size(P,2); % No. mixture components

SRows = zeros(n*nv,ni+2); % The 1st col: Cluster#, 2nd: ID (seq. number)

for v=1:nv
    S = mnrnd(n,P(:,v)');

    [IDX] = find(S > 0);
    NUM = S(IDX);

    pos = (n*(v-1)+1);
    for i=1:length(IDX)
        SRows(pos:(pos + NUM(i) - 1),3:end) = repmat(PermI(IDX(i),:),[NUM(i),1]);
        pos = pos + NUM(i);
    end
    SRows((n*(v-1)+1):(n*v),1) = v;
    SRows((n*(v-1)+1):(n*v),2) = 1:n;
end;

%%
% Write the samples to file
outpath = 'data/MATLAB-data';
csvwrite(strcat(outpath,'/PL-mix.csv'),SRows);

N = 5;
K = 3;
G_orig = [
    0.10,0.56,0.30,0.89,0.02;
    0.50,0.10,0.15,0.05,0.20;
    0.25,0.30,0.10,0.10,0.25]';
G_orig = G_orig*diag(1./sum(G_orig,1));
csvwrite(strcat(outpath,'/PL-mix-gamma-K',num2str(K),'_N',num2str(N),'.csv'),G_orig);

%%
% plackmm uses MM to fit the Plackett-Luce model.
% The a input matrix is nx3, where each row contains
%   Column 1:  individual ID (1 through M)
%   Column 2:  contest ID (1 through N)
%   Column 3:  rank 
V_hat = zeros(nv,ni);
for v=1:nv
    RM = zeros(n*ni,3); % The Rank Matrix
    pos = 1;
    for i=1:n
        RM(pos:(pos+ni-1),1) = SRows(n*(v-1)+i,3:end)';
        RM(pos:(pos+ni-1),2) = SRows(n*(v-1)+i,2);
        RM(pos:(pos+ni-1),3) = 1:ni;
        pos = pos + ni;
    end
    V_hat(v,:) = plackmm(RM);
end
V
V_hat
%%

