function [norm,time,swaps,R,H] = LSBI_Det_Symmetric (A,r,m,R)

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

det_ref = abs(det(A(R,R)));

while flag > 0
    
    flag = 0;
    
    flag_change = 0;
    
    for i = (1:r)
 
        R_aux = R;
        
        for j = (1:m-r)

            R_aux(i) = Rb(j);
   
            value_det = abs(det(A(R_aux,R_aux)));
            
            if value_det > det_ref
                
                flag_change = 1;
                
                i_aux = i;
                
                j_aux = j;
                
                det_ref = value_det;
            end
            
            R_aux(i) = R(i);
            
        end
    end
        if flag_change == 1
            
            swaps = swaps + 1;
            
            el_save = R(i_aux);
            
            R(i_aux) = Rb(j_aux);
            
            Rb(j_aux) = el_save;
            
            flag   = flag + 1;
            
            flag_change = 0;
            
        end
    
end

H = inv(A(R,R));

norm = sum(sum(abs(H)));


time = toc;

end




