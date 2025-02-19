function [normLSLAFI_Det,timeLSLAFI_Det,swaps_LSLAFI_Det,R_LSLAFI_Det,C_LSLAFI_Det] = LSLAFI_Det(A,r,m,n,R,C)

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

for j = (1:n)
    
    n1 = find(C==j);
    
    if isempty(n1)==1
        
        Cb = [Cb;j];
    end
end


Ar = A(R,C);  % Starting block
flag = 1;

%% Local Search

while flag > 0 && swaps <= 1000
    
    flag = 0;
    
    % FOR ROWS
    
    [L1,U1,P1] = lu(Ar');
    
    for i = (1:m-r)
        
        % LU Factoration
        
        b = P1 * A(Rb(i),C)';
        
        y = L1\b;
        
        alfa = U1\y;
        
        % Changing the Ar matrix
        
        [biggest_alfa,local_alfa] = max(abs(alfa));                                       % Find the biggest value that increases the determinant, and his local
        
        if abs(biggest_alfa) > 1.000001
            
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
        
        % Fatoração LU
        
        b = P2 * A(R,Cb(i));
        
        y = L2\b;
        
        alfa = U2\y;
        
        % Changing Ar
        
        [biggest_alfa,local_alfa] = max(abs(alfa));
        
        if abs(biggest_alfa) > 1.000001
            
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

timeLSLAFI_Det = toc;
R_LSLAFI_Det=R;
C_LSLAFI_Det=C;
swaps_LSLAFI_Det=swaps;
normLSLAFI_Det = sum(sum(abs(inv(full(A(R,C))))));




end




