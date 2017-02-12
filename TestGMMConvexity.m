% We test the convexity of PL Mixture model under 
% Generalized Method of Moments formulation

clear all

outpath = 'data/MATLAB-data';
inpath  = 'data/MATLAB-data';

Rankings = csvread(strcat(outpath,'/PL-mix.csv'));

rng('default'); % MATLAB

% Just use random 10 rankings for test
t_indx = randsample(size(Rankings,1),50);
Tr = Rankings(t_indx,3:end);
Pd = {};
for i=1:size(Tr,1)
    Pd{i} = convertToPairwise(Tr(i,:));
end

Ns = 400;
N = size(Tr,2); % Number of items being ranked
D = size(Tr,1); % Number of detectors
K = 2; % Number of mixture components
Gs = {};
Zs = {};
Ps = {};
for i = 1:Ns
    Gs{i} = zeros(N,K);
    for k = 1:K
        Gs{i}(:,k) = drchrnd(ones(1,N),1)';
    end
    for k = 1:K
        Ps{i} = drchrnd(ones(1,K),1)';
    end
    Zs{i} = zeros(K,D);
    for d = 1:D
        Zs{i}(:,d) = drchrnd(ones(1,K),1)';
        %Zs{i}(:,d) = emProbPLMembership( Tr(d,:), Gs{i}, Ps{i}, N, K );
    end
end

%%
%rng('default'); % MATLAB

%alpha = 0.3;
tries = 50000;
Results = zeros(1,tries);
for i = 1:tries
    alpha = unifrnd(0,1);
    Pts = randsample(Ns,2);
    G_x = Gs{Pts(1)};
    Z_x = Zs{Pts(1)};
    P_x = Ps{Pts(1)};
    G_y = Gs{Pts(2)};
    Z_y = Zs{Pts(2)};
    P_y = Ps{Pts(2)};
    %[G,Z,P] = getAffinePoint(Tr, G_x,P_x,G_y,P_y,alpha,D,N,K);
    [G,Z] = getFullAffinePoint(G_x,Z_x,G_y,Z_y,alpha,D,N,K);
    f_x = gmmPLMixObjective(Pd,G_x,Z_x,D,N,K);
    f_y = gmmPLMixObjective(Pd,G_y,Z_y,D,N,K);
    f_a = gmmPLMixObjective(Pd,G,Z,D,N,K);
    Results(i) = (f_a > (alpha*f_x + (1-alpha)*f_y));
    if (Results(i) == 1)
        [f_x f_y f_a (alpha*f_x + (1-alpha)*f_y)]
    end
end
sum(Results) % should get zero

%% Test convexity of pi_{kd}
rng('default'); % MATLAB

Gs = {};
Zs = {};
for i = 1:Ns
    Gs{i} = zeros(N,K);
    for k = 1:K
        Gs{i}(:,k) = drchrnd(ones(1,N),1)';
    end
    Zs{i} = zeros(K,D);
    for d = 1:D
        Zs{i}(:,d) = drchrnd(ones(1,K),1)';
    end
end
%alpha = 0.3;
tries = 500;
Results = zeros(1,tries);
for i = 1:tries
    alpha = unifrnd(0,1);
    Pts = randsample(Ns,2);
    G_x = Gs{Pts(1)};
    Z_x = Zs{Pts(1)};
    G_y = Gs{Pts(2)};
    Z_y = Zs{Pts(2)};
    Z = getAffinePointForCoorDesc(Z_x,Z_y,alpha,D,K);
    f_x = gmmPLMixObjective(Pd,G_x,Z_x,D,N,K); % All objectives evaluated at G_x
    f_y = gmmPLMixObjective(Pd,G_x,Z_y,D,N,K);
    f_a = gmmPLMixObjective(Pd,G_x,Z,D,N,K);
    Results(i) = (f_a > (alpha*f_x + (1-alpha)*f_y));
end
sum(Results) % should get zero

%% Test convexity of sum((|phi_k| / |Z_k|) * sum(z_ki*phi_ki))
%rng('default'); % MATLAB

Ns = 400;
tries = 10000;
totalSum = 0;
for ii=1:1
    Results = zeros(1,tries);
    Zs = {};
    for i = 1:Ns
        for d = 1:D
            Zs{i}(:,d) = drchrnd(ones(1,K),1)';
        end
    end
    Phi = unifrnd(0,5,D,K);
    Zt = {};
    for i=1:tries
        Pts = randsample(Ns,2);
        alpha = unifrnd(0,1);
        Zt{1} = Zs{Pts(1)};
        Zt{2} = Zs{Pts(2)};
        Zt{3} = alpha*Zt{1} + (1-alpha)*Zt{2};
        obj = [0;0;0];
        for j=1:3
            Z = Zt{j};
            for k=1:K
                %  This is neither convex nor concave
                %{
                obj(j) = obj(j) + ...
                    sqrt(Phi(:,k)'*Phi(:,k)) ...
                    * (1/sum(Z(k,:))) ...
                    * Z(k,:)*Phi(:,k);
                %}
                %  This is concave
                %{
                obj(j) = obj(j) + ...
                    sqrt(Phi(:,k)'*Phi(:,k)) ...
                    * (1/sqrt((Z(k,:)*Z(k,:)'))) ...
                    * Z(k,:)*Phi(:,k);
                %}
                % This appears convex
                %%{
                obj(j) = obj(j) + ...
                    sqrt(Phi(:,k)'*Phi(:,k)) * sqrt((Z(k,:)*Z(k,:)')) ...
                    * (1/(sum(Z(k,:)))^2) ...
                    * Z(k,:)*Phi(:,k);
                %%}
                % This is neither convex nor concave
                %{
                for m=1:D
                    for n=1:D
                        obj(j) = obj(j) + Z(k,m)*Z(k,n)*Phi(m,k)*Phi(n,k);
                    end
                end
                obj(j) = obj(j) * (1/(sum(Z(k,:)))^2);
                %}
            end
        end
        Results(i) = (obj(3) > (alpha*obj(1) + (1-alpha)*obj(2)));
    end
    totalSum = totalSum + sum(Results);
    [ii sum(Results) totalSum]
end

%%

A = 0:0.01:1;
K = 2;
Phi = unifrnd(0,5,D,K);
%n = size(A,2);
Z = [A' 1-A'];
%B = repmat(A,n,1);
%Z = [reshape(B,[n*n 1]) repmat(A',n,1)];
F = arrayfun(@(y) surrObj(Z(y,:), Phi), 1:size(Z,1), 'UniformOutput', false);
cell2mat(F)

surrObj(Z(1,:)', Phi)
