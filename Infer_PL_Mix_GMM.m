%% Use CVX to solve for the PL mixture parameters

clear all

rng('default'); % MATLAB

outpath = 'report';
%inpath  = 'data/MATLAB-data';
inpath = 'report';
%inpath = 'benchmark';

% Note: This only applies to K=1 (Single underlying PL distribution)
data = csvread(strcat(inpath,'/PL-mix-K1-N30-D20.csv'));
%data = csvread(strcat(inpath,'/shuttle-ranks.csv'));
%{
G_orig = [
    0.10,0.56,0.30,0.89,0.02;
    0.50,0.10,0.15,0.05,0.20;
    0.25,0.30,0.10,0.10,0.25]';
%}
G_orig = csvread(strcat(inpath,'/PL-mix-gamma-K1_N30.csv'));
G_orig = G_orig*diag(1./sum(G_orig,1));

% Use fee random rankings for test
%t_indx = randsample(size(Rankings,1),100);
%Tr = Rankings(t_indx,3:end);
Tr = data(:,3:end);
Pd = {};
for i=1:size(Tr,1)
    Pd{i} = convertToPairwise(Tr(i,:));
end

N = size(Tr,2); % Number of items being ranked
D = size(Tr,1); % Number of detectors
%K = 1; % Number of mixture components
K = size(G_orig,2);

% Initialize gamma and pi's
Gs = zeros(N,K);
for k = 1:K
    Gs(:,k) = drchrnd(ones(1,N),1)';
end
for k = 1:K
    Ps = drchrnd(ones(1,K),1)';
end
Zs = zeros(K,D);
for d = 1:D
    %Zs(:,d) = emProbPLMembership( Tr(d,:), Gs, Ps, N, K );
    Zs(:,d) = drchrnd(ones(1,K),1)';
end

%% Main iteration

% cd cvx/cvx
% cvx_setup
datestr(clock, 0)

Gs_saved = [];
max_epochs = 1;
S = zeros(N,N);
for epoch=1:max_epochs
    % =============================================
    % Infer the g's keeping the z's constant
    % =============================================
    for k = 1:K
        S(:,:) = 0;
        for d = 1:D
            S(:,:) = S(:,:) + Zs(k,d) * Pd{d};
        end
        S(:,:) = S(:,:) / sum(Zs(k,:));
        % compute Gamma
        % cvx to solve g(N, K)
        cvx_begin
            variable g(N,1) %nonnegative
            minimize (g' * (S' * S) * g)
            %minimize (g' * (S' * S + 0.1*eye(N)) * g)
            subject to
                g >= 0
                %sum(g) == 1
        cvx_end
        % end of cvx to solve g(N, K)
        Gs(:,k) = g;
    end
    if (mod(epoch,1) == 0)
        Gs_saved = [Gs_saved;[ones(N,1)*epoch,Gs]];
        fprintf('epoch=%d\n',epoch);
    end
end

datestr(clock, 0)

Gs = Gs*diag(1./sum(Gs,1));
%%
csvwrite(strcat(outpath,'/PL_Gs_K',num2str(K),'_N',num2str(N),'-inferred-GMM.csv'),Gs_saved);
%csvwrite(strcat(outpath,'/PL_Shuttle_Gs_K',num2str(K),'_N',num2str(N),'-inferred-GMM.csv'),Gs_saved);

%%
[MSE, nepochs] = computeMSE(Gs_saved, G_orig, N, K);
plot(1:nepochs, MSE);
