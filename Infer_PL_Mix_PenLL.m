clear all

outpath = 'report';
%inpath  = 'data/MATLAB-data';
inpath = 'report';

data = csvread(strcat(inpath,'/PL-mix-K1-N30-D20.csv'));
G_orig = csvread(strcat(inpath,'/PL-mix-gamma-K1_N30.csv'));
G_orig = G_orig*diag(1./sum(G_orig,1));

%K = 2;
K = size(G_orig,2);
Tr = data(:,3:end);
Pd = {};
for i=1:size(Tr,1)
    Pd{i} = convertToPairwise(Tr(i,:));
end
D = size(Tr,1);
N = size(Tr,2);


%% Test btls

beta = 0.5;
gamma = 0.2;
f  = @(X) X(1)^2+10*(X(2)^2);
df = @(X) [ 2*X(1) ; 20*X(2) ];
[X, k, Fx] = btls( [1;-1], f, df, beta, gamma )

%%

% Initialize Z with Dirichlet samples
rng('default'); % MATLAB
Zs = drchrnd(ones(1,K),D);

% Initialize mixture proportions
Pi = sum(Zs)/D;

% Initialize Gs's to random
Gs_saved = [];
Gs = drchrnd(ones(1,N),K)';
Gs_old = zeros(N,K);

Items = 1:N;
Ordered_ranks = zeros(D,N);
for d = 1:D
    [~,I] = sort(Tr(d,:));
    Ordered_ranks(d,:) = Items(I);
end

S = zeros(N,N,K);
SS = zeros(N,N,K);
lambda = 1; %1/D; % factor for penalty term
beta = 0.5;
gamma = 0.2;
maxepochs = 50;
options = optimoptions('fmincon','Algorithm','interior-point',...
    'GradObj','on','GradConstr','on');
lb = []; ub = [];
for epoch = 1:maxepochs
  
    % M-Step
    % ------
  
    % MLE of Pi
    Pi = sum(Zs)/D;

    for k = 1:K
        S(:,:,k) = 0;
        for d = 1:D
            S(:,:,k) = S(:,:,k) + Zs(d,k) * Pd{d};
        end
        S(:,:,k) = S(:,:,k) / sum(Zs(:,k));
        SS(:,:,k) = S(:,:,k)' * S(:,:,k);
    end
    
    % Compute the MLE of Gs_kn using gradient descent and BTLS
    %Gs_old(:,:) = Gs(:,:);
    %fn = @(Gs) gammaObjFn(Tr, Gs, Gs_old, Zs, SS, lambda);
    %dfn = @(Gs) gammaObjDFn(Tr, Ordered_ranks, Gs, Gs_old, Zs, SS, lambda);
    %[X, k, Fx] = btls( Gs, fn, dfn, beta, gamma );
    
    G = reshape(Gs,1,[]);
    G_old = G;
    fn = @(X) gammaObjFn(Tr, Ordered_ranks, X, G_old, Zs, SS, lambda);
    [x,fval,exitflag] = fmincon(fn,G,[],[],[],[],lb,ub,@confungrad,options);
    
    Gs = reshape(x, N, K);
    
    % E-Step
    % ------
    for d = 1:D
        for k = 1:K
            tz = Pi(k);
            for i = 1:N
                tz = tz * Gs(Tr(d,i),k) / sum(Gs(Tr(d,i:N),k));
            end
            Zs(d,k) = tz;
        end
    end
    Zs = diag(1./ sum(Zs,2))*Zs; % normalize to 1
    
    if (mod(epoch,1) == 0)
        Gs_saved = [Gs_saved;[ones(N,1)*epoch,Gs]];
        fprintf('epoch=%d\n',epoch);
    end
end

Gs = Gs*diag(1./sum(Gs,1));

%%
csvwrite(strcat(outpath,'/PL_Gs_K',num2str(K),'_N',num2str(N),'-inferred-intpoint.csv'),Gs_saved);

%%
[MSE, nepochs] = computeMSE(Gs_saved, G_orig, N, K);
plot(1:nepochs, MSE);
