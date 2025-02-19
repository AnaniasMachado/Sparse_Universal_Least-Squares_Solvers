function [norm,time,swaps,R,H] = LSFI_Norm_Symmetric (A,r,m,R)

tic;

swaps = 0;

%% Knowing the rows and columns that i'm not using.
%  The current rows = R and the current columns = C

Rb = [];

for i = (1:m)
    
    m1 = find(R==i);
    
    if isempty(m1)==1
        
        Rb = [Rb;i];
    end
end

flag = 1;

%% Local Search

norm_ref = sum(sum(abs(inv(A(R,R)))));

while flag > 0
    
    flag = 0;
    
    for i = (1:r)
        
        R_aux = R;
        
        for j = (1:m-r)
            
            R_aux(i) = Rb(j);
            
            if rank(A(R_aux,R_aux))>=r
            %if abs(det(A(R_aux,R_aux))) > 0
                
                value_norm = sum(sum(abs(inv(A(R_aux,R_aux)))));

                if value_norm < norm_ref

                    norm_ref = value_norm;

                    swaps = swaps + 1;

                    el_save = R(i) ;

                    R(i) = Rb(j);

                    Rb(j) = el_save;

                    flag   = flag + 1;

                    break

                else

                    R_aux(i) = R(i);

                end
                
            else
                 R_aux(i) = R(i);
                 
            end
                
            
        end
        
    end

    H = inv(A(R,R));
    
    norm = sum(sum(abs(H)));
    time = toc;
   
    
end




