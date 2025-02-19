function [norm,time,swaps,C_out] = LSFI_Norm_P3(A,r,n,C)

tic

%% Knowing the rows and columns that i'm not using.
%  The current rows = R and the current columns = C

m = size(A,1);
Cb = [];

for j = (1:n)
    
    n1 = find(C==j);
    
    if isempty(n1)==1
        
        Cb = [Cb;j];
    end
end

A_hat = A(:,C);

H_hat = (A_hat'*A_hat) \ A_hat';
% H_hat = pinv(A_hat);
value_ref = sum(sum(abs(H_hat)));

swaps = 0;

flag = 1;

%% LOCAL SEARCH

%% FOR COLUMNS

while flag>0
    
    flag = 0;
    
    for column = (1:n-r)
        
        %flag_c = 0;
        v = H_hat*A(:,Cb(column));
        
        sum_value = 0;
        
        for i = (1:r)
            
            if abs(v(i))<1e-8
                continue;
            end
            
            abs_row_value = sum(abs(H_hat(i,:)));              % valor da soma da linha onde desejamos trocar a coluna (onde tem o vr)
            
            sum_value = sum_value + abs((1/v(i)))*abs_row_value;            % Row of the Vr that will be divide by vr
            
            for j = (1:r)
                
                if j~=i
                    
                    for k = (1:m)
                        
                        sum_value = sum_value + abs(H_hat(j,k) - (v(j)/v(i))*H_hat(i,k));    % A propria linha - v(j)/v(r) * linha do vr
                        
                    end
                    
                end
            end
            
            if sum_value < value_ref
                
                swaps = swaps + 1;
                
                flag = flag+1;
                
                H_hat(i,:) = (1/v(i))*H_hat(i,:);
                
                for j = (1:r)
                    if j ~= i
                        H_hat(j,:) =  H_hat(j,:) - H_hat(i,:)*v(j);
                    end
                end
                
                Cb_save = Cb(column);
                
                Cb(column) = C(i);
                
                C(i) = Cb_save;
                
                value_ref = sum_value;
                
                break                       %Go back to i = 1;
                
            else
                sum_value = 0;
                
            end
            
        end
    end
end

% A_hat = A(:,C);
% H_hat = inv(A_hat'*A_hat) * A_hat';
time = toc;
C_out = C;
norm = value_ref;

end