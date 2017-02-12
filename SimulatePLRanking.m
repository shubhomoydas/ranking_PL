
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
V = [0.1,0.56,0.3,0.89,0.02];
V = V / sum(V);
I = 5:-1:1;
PermI = perms(I); % P(1:10,:)
np = size(PermI,1);
P = zeros(np,1);
for i=1:np
    tpi = PermI(i,:);
    %tpi = [2,5,3,4,1];
    P(i,1) = computePLProb(tpi,V);
end;
%%

rng('default'); % MATLAB
n = 100;
ni = size(PermI,2);
S = mnrnd(n,P');

[IDX] = find(S > 0);
NUM = S(IDX);

SRows = zeros(n,ni+1); % The first col is ID -- a sequential number
pos = 1;
for i=1:length(IDX)
    SRows(pos:(pos + NUM(i) - 1),2:end) = repmat(PermI(IDX(i),:),[NUM(i),1]);
    pos = pos + NUM(i);
end
SRows(:,1) = 1:n;

% plackmm uses MM to fit the Plackett-Luce model.
% The a input matrix is nx3, where each row contains
%   Column 1:  individual ID (1 through M)
%   Column 2:  contest ID (1 through N)
%   Column 3:  rank 
RM = zeros(n*ni,3); % The Rank Matrix
pos = 1;
for i=1:n
    RM(pos:(pos+ni-1),1) = SRows(i,2:end)';
    RM(pos:(pos+ni-1),2) = SRows(i,1);
    RM(pos:(pos+ni-1),3) = 1:ni;
    pos = pos + ni;
end
%%

v_hat = plackmm(RM);
