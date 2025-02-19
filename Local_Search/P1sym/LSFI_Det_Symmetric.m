function [norm,time,swaps,R,H] = LSFI_Det_Symmetric (A,r,m,R)

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
    
    for i = (1:r)
        
        R_aux = R;

        for j = (1:m-r)
            
            R_aux(i) = Rb(j);
            
            value_det = abs(det(A(R_aux,R_aux)));
            
            if value_det > det_ref
                
                det_ref = value_det;
                
                swaps = swaps + 1;
                
                el_save = R(i) ;
                
                R(i) = Rb(j);
                
                Rb(j) = el_save;
                
                flag   = flag + 1;
                
                break
                
            else
                
                R_aux(i) = R(i);
                
            end
            
        end
        
    end

    H = inv(A(R,R));
    
    norm = sum(sum(abs(H)));
    time = toc;
   
    
end




