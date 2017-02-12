function [ MSE, nepochs ] = computeMSE( Gs_saved, G_orig, N, K )
%COMPUTEMSE Summary of this function goes here
%   Detailed explanation goes here

nepochs = max(Gs_saved(:,1));
TGs = {};
MSE = zeros(nepochs,1);
for epoch=1:nepochs
    Gs = Gs_saved(Gs_saved(:,1)==epoch,2:end);
    % ============================================
    % Find the pair-wise similar clusters between
    % inferred and original
    % ============================================
    G = zeros(N,K);
    Sim = zeros(K,K);
    GsNorm = Gs*diag(1./sum(Gs,1));
    for k=1:K
        for kk=1:K
            Sim(k,kk) = sum((GsNorm(:,k)-G_orig(:,kk)).^2);
        end
    end
    % Step 1: Find two closest matching pairs from the
    %         all-pairs distance matrix Sim
    % Step 2: Pair up the above two clusters and remove
    %         from further consideration by removing
    %         them from corresponding arrays (Kx, Ky, Sim)
    % Step 3: If pairs still remain to be matched in Kx/Ky, go to Step 1
    Kx = 1:K;
    Ky = 1:K;
    for k=1:K
        smin = Inf;
        id_x = 0;
        id_y = 0;
        ksim = size(Sim,1);
        for i=1:ksim
            [v,id] = min(Sim(i,:));
            if (v < smin)
                smin = v;
                id_x = id;
                id_y = i;
            end;
        end
        G(:,Kx(id_x)) = GsNorm(:,Ky(id_y));
        MSE(epoch) = MSE(epoch) + sum((G(:,Kx(id_x))-G_orig(:,Kx(id_x))).^2);
        Kx(id_x) = [];
        Ky(id_y) = [];
        Sim(:,id_x) = [];
        Sim(id_y,:) = [];
    end
    % ====

    TGs{epoch} = G;
    
end

%plot(1:nepochs,MSE);

end

