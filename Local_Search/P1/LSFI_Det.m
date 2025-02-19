function [normLSFI_Det,timeLSFI_Det,swaps_LSFI_Det,R_LSFI_Det,C_LSFI_Det] = LSFI_Det (A,r,m,n,R,C)

tic;

swaps = 0;

%% Knowing the rows and columns that i'm not using. 
%  The current rows = R and the current columns = C

Rb = [];                                                                            
Cb = [];                                                                            

for i = (1:m)
    
    m1 = find(R==i);
    
    if isempty(m1)==1
        
        Rb = [Rb;i];
    end
end

Ar = A(R,C);

for j = (1:n)
    
    n1 = find(C==j);
    
    if isempty(n1)==1
        
        Cb = [Cb;j];
    end
end


flag = 1;

%% Local Search

while flag > 0  && swaps <= 1000
    
    flag = 0;
    
    % FOR ROWS
    
    [L1,U1,P1] = lu(Ar');
    
    for i = (1:m-r)
        
        % LU Factoration
        
        b = P1 * A(Rb(i),C)';
        
        y = L1\b;
        
        alfa = U1\y;
        
        % Changing Ar
        
        local_alfa = find(abs(alfa)>1.000001);
        
        T = isempty(local_alfa);
        
        if T==0
            
            local_alfa = local_alfa(1);
            
            swaps = swaps + 1;
            
            el_save = R(local_alfa) ;
            
            R(local_alfa) = Rb(i);
            
            Rb(i) = el_save;
            
            flag   = flag + 1;
            
            Ar = A(R,C);
            
            [L1,U1,P1] = lu(Ar');
        end
    end
    
    % FOR COLUMNS
    
    [L2,U2,P2] = lu(Ar);
    
    for i = (1: n-r)
        
        % LU Factoration
        
        b = P2 * A(R,Cb(i));
        
        y = L2\b;
        
        alfa = U2\y;
        
        
        % Changing Ar
        
        local_alfa = find(abs(alfa)>1.000001);
        
        T = isempty(local_alfa);
        
        if T==0
            
            local_alfa = local_alfa(1);
            
            swaps = swaps + 1;
            
            el_save = C(local_alfa) ;
            
            C(local_alfa) = Cb(i);
            
            Cb(i) = el_save;
            
            flag   = flag + 1;
            
            Ar = A(R,C);
            
            [L2,U2,P2] = lu(Ar);
        end
        
    end
    
%     if det(Ar) == 0
%         flag = flag+1;
%     end
    
end

timeLSFI_Det = toc;
R_LSFI_Det=R;
C_LSFI_Det=C;
swaps_LSFI_Det=swaps;
normLSFI_Det = sum(sum(abs(inv(A(R,C)))));


end