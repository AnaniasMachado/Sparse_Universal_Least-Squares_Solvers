function [norm,time,swaps,C_out] = LSBI_Det_P3 (A,r,n,R,C)

tic;

swaps = 0;

%% Knowing the rows and columns that i'm not using.
%  The current rows = R and the current columns = C


Cb = [];


for j = (1:n)
    
    n1 = find(C==j);
    
    if isempty(n1)==1
        
        Cb = [Cb;j];
    end
end


% Ar = A(R,C);  % Starting block

flag = 1;
%% Local Search

while flag> 0
    
    Ar = A(R,C);
    
    flag = 0;
    
    % FOR COLUMNS
    
    [L2,U2,P2] = lu(Ar);
    biggest_alfa = 1 ;
    
    for i = (1: n-r)
        
        % Fatoração LU
        
        b = P2 * A(R,Cb(i));
        
        y = L2\b;
        
        alfa = U2\y;
        
        % Changing Ar
        [new_biggest_alfa,new_local_alfa] = max(abs(alfa));
        
        if (new_biggest_alfa > biggest_alfa)
            biggest_alfa = new_biggest_alfa ;
            local_alfa = new_local_alfa ;
            swap_column = i;
        end
    end
    if biggest_alfa > 1
            
        swaps = swaps + 1;

        el_save = C(local_alfa) ;

        C(local_alfa) = Cb(swap_column);

        Cb(swap_column) = el_save;

%         Ar = A(R,C);

%         [L2,U2,P2] = lu(Ar);

        flag = flag + 1;
    end
        
    
end

% End of cicle for
time = toc;
A_hat =  A(:,C);

H_hat = (A_hat'*A_hat) \ A_hat';
% H_hat = pinv(A_hat);


norm = sum(sum(abs(H_hat)));
C_out = C;


end


