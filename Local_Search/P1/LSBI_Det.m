function [normLSBI_Det,timeLSBI_Det,swaps_LSBI_Det,R_LSBI_Det,C_LSBI_Det] = LSBI_Det (A,r,m,n,R,C)

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

while flag > 0  && swaps <= 1000
    
    flag = 0;
    
    % FOR ROWS
    
    [L1,U1,P1] = lu(Ar');
    biggest_alfa_row = 1 ;
    
    for i = (1:m-r)
        
        % LU Factoration
        
        b = P1 * A(Rb(i),C)';
        
        y = L1\b;
        
        alfa = U1\y;
        
        % Changing the Ar matrix
        
        [new_biggest_alfa,new_local_alfa] = max(abs(alfa));                                       % Find the biggest value that increases the determinant, and his local
        
        if new_biggest_alfa> biggest_alfa_row
            biggest_alfa_row = new_biggest_alfa ;
            local_alfa_row = new_local_alfa ;
            swap_row = i;
        end
    end
    
    % FOR COLUMNS
    
    [L2,U2,P2] = lu(Ar);
    biggest_alfa_col = 1 ;
    
    for i = (1: n-r)
        
        % Fatoração LU
        
        b = P2 * A(R,Cb(i));
        
        y = L2\b;
        
        alfa = U2\y;
        
        % Changing Ar
        
        [new_biggest_alfa,new_local_alfa] = max(abs(alfa));
        
       if (new_biggest_alfa > biggest_alfa_col)
            biggest_alfa_col = new_biggest_alfa ;
            local_alfa_col = new_local_alfa ;
            swap_column = i;
        end
    end
    if max(biggest_alfa_row,biggest_alfa_col) > 1.000001
        
    if biggest_alfa_row > biggest_alfa_col
            
            swaps = swaps + 1;
            
            el_save = R(local_alfa_row) ;
            
            R(local_alfa_row) = Rb(swap_row);
            
            Rb(swap_row) = el_save;
            
            flag   = flag + 1;
            
            Ar = A(R,C);
            
    %        [L1,U1,P1] = lu(Ar');
%     end
%     
%     if biggest_alfa_col > 1
    else
            
            swaps = swaps + 1;
            
            el_save = C(local_alfa_col) ;
            
            C(local_alfa_col) = Cb(swap_column);
            
            Cb(swap_column) = el_save;
            
            flag   = flag + 1;
            
            Ar = A(R,C);
            
            [L2,U2,P2] = lu(Ar);
    end
    end
        
    
%     if det(Ar) == 0
%         flag = flag+1;
%     end
%     
end


%log_det = log(abs(det(Ar)));

% H = inv(Ar);
% norm = sum(sum(abs(H)))
% nswaps = swaps
% time = toc
% 
timeLSBI_Det = toc;
R_LSBI_Det=R;
C_LSBI_Det=C;
swaps_LSBI_Det=swaps;
normLSBI_Det = sum(sum(abs(inv(A(R,C)))));
end




