function [ DGs ] = gammaObjDFn( Tr, Ordered_ranks, Gs, Gs_old, Zs, SS, lambda )
%GAMMAOBJDFN Summary of this function goes here
%   Detailed explanation goes here
    
    D = size(Tr,1);
    N = size(Tr,2);
    K = size(Zs,2);
    
    DGs = zeros(1,N*K);
    for k=1:K
        Gs_k = Gs((N*(k-1)+1):(N*k));
        Gs_old_k = Gs_old((N*(k-1)+1):(N*k));
        idx = (k-1)*N;
        for n=1:N
            % we assume that we get complete rankings on N items from
            % each detector. Therefore, we can simplify the indicator
            % function in the numerator of the derivation.
            nv = sum(Zs(:,k));
            tv = 0;
            for d=1:D
                s_l = 0;
                for i = 1:N
                    %[d n size(Ordered_ranks)]
                    if (Ordered_ranks(d,n) >= i) % this check is equivalent to the 'which'
                        s_l = s_l + (1/sum(Gs_old_k(Tr(d,i:N))));
                    end
                end
                tv = tv + Zs(d,k)*s_l;
            end
            %DGs(n,k) = -((nv/Gs(n,k)) - tv - 100/Gs(n,k)); % maximize
            DGs(idx+n) = -((nv/Gs_k(n)) - tv); % maximize
        end
        % penalty term derivative
        DGs((idx+1):(idx+N)) = DGs((idx+1):(idx+N)) + (lambda*SS(:,:,k)*Gs_k')';
    end
    
end

