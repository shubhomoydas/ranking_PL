
clear all

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

%outpath = 'data/MATLAB-data';
outpath = 'report';

N = 30; % No. of items
K =  3; % No. of clusters

Pi = ones(1,K) / K; % Distribution of clusters
D  = 20; % No. of competitions

nD = ceil(Pi*D);
nD(K) = D - sum(nD(1:(K-1))); % make sure the sum = D

%% Generate the Gamma parameters

rng('default'); % MATLAB
Gamma = drchrnd(ones(1,N),K)'; % Gamma parameters
% sum(Gamma) ==> 1
csvwrite(strcat(outpath,'/PL-mix-gamma-K',num2str(K),'_N',num2str(N),'.csv'),Gamma);

%%

rng('default'); % MATLAB

SRows = zeros(D,N+2); % The 1st col: Cluster#, 2nd: ID (seq. number)
d = 1;
for k=1:K
    for i=1:nD(k)
        SRows(d,3:end) = sampleFromPLDistribution(Gamma(:,k));
        SRows(d,1:2) = [k d];
        d = d+1;
    end
    fprintf('Samples Generated for K=%d\n',k);
end

%%
% Write the samples to file
csvwrite(strcat(outpath,'/PL-mix-K',num2str(K),'-N',num2str(N),'-D', ...
    num2str(D),'.csv'),SRows);

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

