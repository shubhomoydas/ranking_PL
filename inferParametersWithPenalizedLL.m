%% Use CVX to solve for the PL mixture parameters

clear all

rng('default'); % MATLAB

outpath = 'data/MATLAB-data';
inpath  = 'data/MATLAB-data';

%Rankings = csvread(strcat(outpath,'/PL-mix-K3-N50-D600.csv'));
%Rankings = csvread(strcat(outpath,'/PL-mix-K1-N50-D600.csv'));
Rankings = csvread(strcat(outpath,'/PL-mix.csv'));

%%{
G_orig = [
    0.10,0.56,0.30,0.89,0.02;
    0.50,0.10,0.15,0.05,0.20;
    0.25,0.30,0.10,0.10,0.25]';
%%}
%G_orig = csvread(strcat(outpath,'/PL-mix-gamma-K1.csv'));
G_orig = G_orig*diag(1./sum(G_orig,1));

% Just use random 100 rankings for test
%t_indx = randsample(size(Rankings,1),100);
%Tr = Rankings(t_indx,3:end);
Tr = Rankings(:,3:end);
Pd = {};
for i=1:size(Tr,1)
    Pd{i} = convertToPairwise(Tr(i,:));
end

N = size(Tr,2); % Number of items being ranked
D = size(Tr,1); % Number of detectors
K = 3; % Number of mixture components

% Initialize gamma and pi's
Gs = zeros(N,K);
for k = 1:K
    Gs(:,k) = drchrnd(ones(1,N),1)';
end
Ps = drchrnd(ones(1,K),1)';
Zs = zeros(K,D);
for d = 1:D
    %Zs(:,d) = emProbPLMembership( Tr(d,:), Gs, Ps, N, K );
    Zs(:,d) = drchrnd(ones(1,K),1)';
end

%% Main iteration

% cd /nfs/guille/u2/d/dassh/classes/CS556/code/cvx/cvx
% cvx_setup

max_epochs = 20;
S = zeros(N,N,K);
SS = zeros(N,N,K);
lambda = 0.1; %1/D; % factor for penalty term
Gs_saved = [];
P_saved = [];
for epoch=1:max_epochs
    
    % =============================================
    % Infer the p's and g's keeping the z's constant
    % =============================================
    Ps = sum(Zs,2)/D;
    for k = 1:K
        S(:,:,k) = 0;
        for d = 1:D
            S(:,:,k) = S(:,:,k) + Zs(k,d) * Pd{d};
        end
        S(:,:,k) = S(:,:,k) / sum(Zs(k,:));
        SS(:,:,k) = S(:,:,k)' * S(:,:,k);
    end
    % compute Gamma
    % cvx to solve g(N, K)
    cvx_begin
        variables g(N,K) %nonnegative
        expressions obj pen(K)
        % first the penalizing term
        for k=1:K
            pen(k) = g(:,k)' * SS(:,:,k) * g(:,k);
        end
        for d=1:D
            for k=1:K
                %obj = obj + Zs(k,d)*p(k);
                for i=1:N
                    obj = obj + Zs(k,d)*log(g(Tr(d,i),k));
                    obj = obj - Zs(k,d)*(sum(g(Tr(d,i:N),k))/sum(Gs(Tr(d,i:N),k)));
                end
            end
        end
        maximize obj - lambda*sum(pen)
        subject to
            g >= 0
            %p >= 0
            sum(g) == 1
            %sum(p) == 1
    cvx_end
    % end of cvx to solve g(N, K)
    
    Gs(:,:) = g;
    Gs_saved = [Gs_saved;[ones(N,1)*epoch,Gs]];
    P_saved = [P_saved;[epoch,Ps']];

    % =============================================
    % Infer the z's keeping the g's constant
    % Since the problem is not convex, we are now
    % going to solve this using a surrogate convex
    % function and interior point method.
    % =============================================

    % The below E-Step for Zs
    for d = 1:D
        Zs(:,d) = emProbPLMembership( Tr(d,:), Gs, Ps, N, K );
    end
    
    fprintf('Completed Epoch %d\n',epoch);
    
end

%%
% Check how the regularizers fared
for k=1:K
    Gs(:,k)'*SS(:,:,k)*Gs(:,k)
end

%%
csvwrite(strcat(outpath,'/PL_Gs_K',num2str(K),'-inferred.csv'),Gs_saved);
csvwrite(strcat(outpath,'/PL_Ps_K',num2str(K),'-inferred.csv'),P_saved);

%% Analyze Convergence
