function [normLSFI_Norm,timeLSFI_Norm,swaps_LSFI_Norm,R_LSFI_Norm,C_LSFI_Norm] = LSFI_Norm(A,r,m,n,R,C)

tic

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

flag = 1;

Br_inv = inv(A(R,C));

swaps = 0;

%% LOCAL SEARCH

while flag > 0  && swaps <= 1000
    
    flag = 0;
    
    %% FOR COLUMNS
    
    for column = (1:n-r)
        
        v = Br_inv*A(R,Cb(column));                  % v  = inv(Br)*column vector we wanna change
        
        abs_value = sum(sum(abs(Br_inv)));           % absolute value of the norm of the inverse
        
        sum_value = 0;
        
        % OBS: Dont forget, you can't do the sum directly of all the rows
        % less the vr, because is in abs, the when you are subtracting,
        % could increase the norm if the value is negative already. And we
        % don't want this.
        
        for i = (1:r)
            
            abs_row_value = sum(abs(Br_inv(i,:)));              % valor da soma da linha onde desejamos trocar a coluna (onde tem o vr)
            
            sum_value = sum_value + abs((1/v(i)))*abs_row_value;            % Row of the Vr that will be divide by vr
            
            for j = (1:r)
                
                if j~=i
                    
                    for k = (1:r)
                        
                        sum_value = sum_value + abs(Br_inv(j,k) - (v(j)/v(i))*Br_inv(i,k));    % A propria linha - v(j)/v(r) * linha do vr
                        
                    end
                    
                end
            end
            
            if sum_value < abs_value - 0.000001
                
                swaps = swaps + 1;
                
                flag = flag+1;
                
                Br_inv(i,:) = (1/v(i))*Br_inv(i,:);
                
                for j = (1:r)
                    if j ~= i
                        Br_inv(j,:) =  Br_inv(j,:) - Br_inv(i,:)*v(j);
                    end
                end
                
                Cb_save = Cb(column);
                
                Cb(column) = C(i);
                
                C(i) = Cb_save;
                
                abs_value = sum_value;
                
                break                       %Go back to i = 1;
                
            else
                sum_value = 0;
                
            end
            
        end
    end
    
    %% FOR ROWS
    
    A = A';
    
    Br_inv = Br_inv';
    
    for row = (1:m-r)

        v = Br_inv*A(C,Rb(row));                      % v  = inv(Br)*vetor coluna que queremos trocar
        
        abs_value = sum(sum(abs(Br_inv)));           %valor da soma em modulo da inversa
        
        sum_value = 0;
        
        for i = (1:r)
            
            abs_row_value = sum(abs(Br_inv(i,:)));              % valor da soma da linha onde desejamos trocar a coluna (onde tem o vr)
            
            sum_value = sum_value + abs((1/v(i)))*abs_row_value;            % Linha do Vr vai ser ela mesma dividido por 1/vr
            
            for j = (1:r)
                
                if j~=i
                    
                    for k = (1:r)
                        
                        sum_value = sum_value + abs(Br_inv(j,k) - (v(j)/v(i))*Br_inv(i,k));    % A propria linha - v(j)/v(r) * linha do vr
                    end
                    
                end
            end
            
            if sum_value < abs_value  -  0.000001
                
                flag = flag+1;
                
                swaps = swaps + 1;
                
                Br_inv(i,:) = (1/v(i))*Br_inv(i,:);
                
                for j = (1:r)
                    if j ~= i
                        Br_inv(j,:) =  Br_inv(j,:) - Br_inv(i,:)*v(j);
                    end
                end
                
                Lb_save = Rb(row);
                
                Rb(row) = R(i);
                
                R(i) = Lb_save;
                
                abs_value = sum_value;
                
                break                       %Go back to i = 1;
                
            else
                sum_value = 0;
                
            end
            
        end
    end
    
    A = A';
    
    Br_inv = Br_inv';
    
end

%H = Br_inv;

%norm = abs_value;


timeLSFI_Norm = toc;
R_LSFI_Norm=R;
C_LSFI_Norm=C;
swaps_LSFI_Norm=swaps;
normLSFI_Norm = abs_value; %sum(sum(abs(inv(full(A(R,C))))));

end