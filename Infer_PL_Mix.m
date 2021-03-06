clear all

%outpath = 'report';
%inpath  = 'data/MATLAB-data';
%inpath = 'report';

%{
outpath = 'report';
inpath = 'benchmark';
data = csvread(strcat(inpath,'/shuttle-ranks.csv'));
K = 2;
%}

%%{
outpath = 'results/PL-ranking';
inpath = '.';
data = csvread(strcat(inpath,'/PL-mix.csv'));
K = 3;
%%}

%{
data = csvread(strcat(inpath,'/PL-mix-K3-N30-D20.csv'));
G_orig = csvread(strcat(inpath,'/PL-mix-gamma-K3_N30.csv'));
G_orig = G_orig*diag(1./sum(G_orig,1));
K = size(G_orig,2);
%}

S = data(:,3:size(data,2));
D = size(S,1);
N = size(S,2);

%%

% Initialize Z with Dirichlet samples
rng('default'); % MATLAB

datestr(clock, 0)
fprintf('Started inference for PL dist without prior, N=%d, K=%d, D=%d\n', N, K, D);

%Z = drchrnd(ones(1,K),D);
Z = ones(D,K)*(1/K);

% Initialize mixture proportions
Pi = sum(Z)/D;

% Initialize Gs's to random
Gs_saved = [];
%Gs = drchrnd(ones(1,N),K)';
Gs = ones(N,K)*(1/N);
Gs_tmp = zeros(N,K);

Items = 1:N;
Ordered_ranks = zeros(D,N);
for d = 1:D
    [~,I] = sort(S(d,:));
    Ordered_ranks(d,:) = Items(I);
end

Gs_sum = zeros(D,N);
ll_prev = -realmax('single');
ll_diff = ll_prev;
maxepochs = 50;
LL = [];
for epoch = 1:maxepochs
    
    % M-Step
    % ------
    
    % MLE of Pi
    Pi = sum(Z)/D;
    
    % MLE of v_kn
    for k = 1:K
        for d = 1:D
            Gs_sum(d,:) = cumsum(Gs(S(d,N:-1:1),k));
            Gs_sum(d,:) = Gs_sum(d,N:-1:1);
        end
        for n = 1:N
            % we assume that we get complete rankings on N items from
            % each detector. Therefore, we can simplify the indicator
            % function in the numerator of the derivation.
            nv = sum(Z(:,k));
            %for d = 1:D {
            %  tv <- tv + Z[d,k]*sum(<Indicator function>)
            %}
            tv = 0;
            for d = 1:D
                s_l = 0;
                for i = 1:N
                    if (Ordered_ranks(d,n) >= i) % this check is equivalent to the 'which'
                        %s_l = s_l + (1/sum(Gs(S(d,i:N),k)));
                        s_l = s_l + (1/Gs_sum(d,i));
                    end
                end
                tv = tv + Z(d,k)*s_l;
            end
            Gs_tmp(n,k) = nv / tv;
        end
    end
    Gs(:,:) = Gs_tmp(:,:);
    
    % E-Step
    % ------
    if (K > 1)
        % TODO: We need to use logarithms for product and normalize for large
        % N
        for k = 1:K
            for d = 1:D
                Gs_sum(d,:) = cumsum(Gs(S(d,N:-1:1),k));
                Gs_sum(d,:) = Gs_sum(d,N:-1:1);
                %tz = Pi(k);
                tz = log(Pi(k));
                for i = 1:N
                    %%tz = tz * Gs(S(d,i),k) / sum(Gs(S(d,i:N),k));
                    %tz = tz + log(Gs(S(d,i),k) / sum(Gs(S(d,i:N),k)));
                    tz = tz + log(Gs(S(d,i),k) / Gs_sum(d,i));
                end
                Z(d,k) = tz;
            end
        end
        %Z = diag(1./ sum(Z,2))*Z; % normalize to 1
        for d = 1:D
            Z(d,:) = Z(d,:) - mean(Z(d,:));
            mx = max(Z(d,:));
            if (mx > 700)
                Z(d,:) = Z(d,:) - (mx - 700);
            end
            Z(d,:) = exp(Z(d,:));
            Z(d,:) = Z(d,:) / sum(Z(d,:)); % normalize to 1
        end
    end
    
    % compute log-likelihood for stopping criteria
    ll = computePLLogLik(S, Gs, Pi, Z);
    LL = [LL; ll];
    ll_diff = ll-ll_prev;
    
    if (mod(epoch,1) == 0)
        Gs_saved = [Gs_saved;[ones(N,1)*epoch,Gs]];
        fprintf('epoch=%d, ll_prev = %f, ll = %f, ll diff = %f\n',...
            epoch, ll_prev, ll, ll_diff);
    end
    
    ll_prev = ll;
    if abs(ll_diff) < 1e-2
        fprintf('Exiting at epoch=%d, ll_prev = %f, ll = %f, ll diff = %f\n',...
            epoch, ll_prev, ll, ll_diff);
        break;
    end
    
end

datestr(clock, 0)
fprintf('Completed inference for PL dist without prior, N=%d, K=%d, D=%d\n', N, K, D);
%Gs = Gs*diag(1./sum(Gs,1));

%{
%csvwrite(strcat(outpath,'/PL_Gs_K',num2str(K),'_N',num2str(N),'-inferred-noprior.csv'),Gs_saved);
%csvwrite(strcat(outpath,'/PL_Gs_K',num2str(K),'_N',num2str(N),'-inferred-noprior-LL.csv'),LL);
csvwrite(strcat(outpath,'/PL_Shuttle_Gs_K',num2str(K),'_N',num2str(N),'-inferred-noprior.csv'),Gs_saved);
csvwrite(strcat(outpath,'/PL_Shuttle_Gs_K',num2str(K),'_N',num2str(N),'-inferred-noprior-LL.csv'),LL);
%}
%%{
csvwrite(strcat(outpath,'/PL_Mix_Gs_K',num2str(K),'_N',num2str(N),'-inferred-noprior.csv'),Gs_saved);
csvwrite(strcat(outpath,'/PL_Mix_Gs_K',num2str(K),'_N',num2str(N),'-inferred-noprior-LL.csv'),LL);
%%}
fprintf('Written to File: Params for PL dist without prior, N=%d, K=%d, D=%d\n', N, K, D);
%%
[MSE, nepochs] = computeMSE(Gs_saved, G_orig, N, K);
plot(1:nepochs, MSE);

