function [ obj, ObjGrad ] = gammaObjFn( Tr, Ordered_ranks, Gs, Gs_old, Zs, SS, lambda )
%GAMMAOBJFN Summary of this function goes here
%   Detailed explanation goes here
    
    D = size(Tr,1);
    N = size(Tr,2);
    K = size(Zs,2);
    
    obj = 0;
    for k=1:K
        Gs_k = Gs((N*(k-1)+1):(N*k));
        Gs_old_k = Gs_old((N*(k-1)+1):(N*k));
        for d=1:D
            for i=1:N
                obj = obj + Zs(d,k)*log(Gs_k(Tr(d,i)));
                obj = obj - Zs(d,k)*(sum(Gs_k(Tr(d,i:N)))/sum(Gs_old_k(Tr(d,i:N))));
            end
        end
        % penalty term
        obj = obj + 0.5*lambda * (Gs_k * SS(:,:,k) * Gs_k');
    end
    
    obj = -obj; % maximize
    
    if nargout > 1
        %fprintf('Returning grad...\n');
        ObjGrad = gammaObjDFn( Tr, Ordered_ranks, Gs, Gs_old, Zs, SS, lambda );
    end
    
end

