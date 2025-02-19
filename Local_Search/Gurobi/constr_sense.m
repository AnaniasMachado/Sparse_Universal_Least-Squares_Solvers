% constr_sense.m

% NOTE: Whether a given constraint set is included or not, that is taken
%   care of by providing the number of rows of the product tightening of P1
%   (nr_PT) and nonsymmetric P2 (nr_NS), as a value of zero (0) for either of
%   these denotes an omission of that type of constraint.

function[sense_form] = constr_sense(m,n,A,nr,constr_type,set_prop)
%keyboard

% Initialize constraint form:
sense_form = cell(1,nr);

count = 0; % Keep track of # of constraint senses stored

% 1-norm constraint formulation: -Z <= H <= Z

for k=1:n*m
    
    count = count + 1;
    
    % H + Z >= 0
    sense_form{1,count} = '>';
    
    count = count + 1;
    
    % Z - H >= 0
    sense_form{1,count} = '>';
    
end % end of for loop, for k=1:n*m

% P1 : AHA - A = 0;
if(strcmp(set_prop{1,1},'included'))
    
    % SVD reduced P1 representation: (SV')H(US) = U'AV (r x r constraints)
    
    r = rank(A);
    
    nr_P1svd = r*r;
    
    for k=1:nr_P1svd
        
        count = count + 1;
    
        if(strcmp(constr_type{1,1},'equality'))
        
            sense_form{1,count} = '=';
            
        elseif(strcmp(constr_type{1,1},'lgr_inequality'))
            
            sense_form{1,count} = '>';
            
        elseif(strcmp(constr_type{1,1},'rgl_inequality'))
            
            sense_form{1,count} = '<';
           
        end % end of if/elseif statements
        
    end % end of for loop, for k=1:nr_P1svd
    
%     % Full P1 representation: AHA - A == 0
%     
%     nr_P1s = m*n;
%     
%     for k=1:nr_P1s
%         
%         count = count + 1;
%         
%         if(strcmp(constr_type{1,1},'equality'))
%         
%             sense_form{1,count} = '=';
%             
%         elseif(strcmp(constr_type{1,1},'lgr_inequality'))
%             
%             sense_form{1,count} = '>';
%             
%         elseif(strcmp(constr_type{1,1},'rgl_inequality'))
%             
%             sense_form{1,count} = '<';
%            
%         end % end of if/elseif statements
%         
%     end
    
else
    
    nr_P1svd = 0;
    
%    nr_P1s = 0;
end

% P3 : AH - (AH)' = 0
if(strcmp(set_prop{1,2},'included'))
    
    nr_P3 = (m^2 - m)/2;
    
    for k=1:nr_P3
        
        count = count + 1;
        
        if(strcmp(constr_type{1,2},'equality'))
        
            sense_form{1,count} = '=';
            
        elseif(strcmp(constr_type{1,2},'lgr_inequality'))
            
            sense_form{1,count} = '>';
            
        elseif(strcmp(constr_type{1,2},'rgl_inequality'))
            
            sense_form{1,count} = '<';
           
        end % end of if/elseif statements
        
    end
    
else
    nr_P3 = 0;
end

% P4 : HA - (HA)' = 0
if(strcmp(set_prop{1,3},'included'))
    
    nr_P4 = (n^2 - n)/2;
    
    for k=1:nr_P4
        
        count = count + 1;
        
        if(strcmp(constr_type{1,3},'equality'))
        
            sense_form{1,count} = '=';
            
        elseif(strcmp(constr_type{1,3},'lgr_inequality'))
            
            sense_form{1,count} = '>';
            
        elseif(strcmp(constr_type{1,3},'rgl_inequality'))
            
            sense_form{1,count} = '<';
           
        end % end of if/elseif statements
    
    end
    
else
    nr_P4 = 0;
end

% sym : H - (H)' = 0
if(strcmp(set_prop{1,10},'included'))
    
%%%%%    nr_P_sym = (m^2 - m)/2;
    nr_P_sym = (m^2 - m);
    
    for k=1:nr_P_sym
        
        count = count + 1;
        
        if(strcmp(constr_type{1,8},'equality'))
        
            sense_form{1,count} = '=';
            
        elseif(strcmp(constr_type{1,8},'lgr_inequality'))
            
            sense_form{1,count} = '>';
            
        elseif(strcmp(constr_type{1,8},'rgl_inequality'))
            
            sense_form{1,count} = '<';
           
        end % end of if/elseif statements
        
    end
    
else
    nr_P_sym = 0;
end


% P2L : HAA+ - H = 0
if(strcmp(set_prop{1,11},'included'))
    
%%%%%    nr_P2L = n*m;
    nr_P2L = m*n;
    
    for k=1:nr_P2L
        
        count = count + 1;
        
        if(strcmp(constr_type{1,9},'equality'))
        
            sense_form{1,count} = '=';
            
        elseif(strcmp(constr_type{1,9},'lgr_inequality'))
            
            sense_form{1,count} = '>';
            
        elseif(strcmp(constr_type{1,9},'rgl_inequality'))
            
            sense_form{1,count} = '<';
           
        end % end of if/elseif statements
        
    end
    
else
    nr_P2L = 0;
end



% Nonsymmetric P2 : HAH - H = 0
if(strcmp(set_prop{1,4},'included'))
    
    nr_NS = n*m;
    
    for k=1:nr_NS
        
        count = count + 1;
        
        if(strcmp(constr_type{1,4},'equality'))
        
            sense_form{1,count} = '=';
            
        elseif(strcmp(constr_type{1,4},'lgr_inequality'))
            
            sense_form{1,count} = '>';
            
        elseif(strcmp(constr_type{1,4},'rgl_inequality'))
            
            sense_form{1,count} = '<';
           
        end % end of if/elseif statements
        
    end
    
else
    nr_NS = 0;
end

% Product Tightening of P1
if(strcmp(set_prop{1,5},'included'))
    
    nr_PT1 = n*m;
    
    for k=1:nr_PT1
        
        count = count + 1;
        
        if(strcmp(constr_type{1,5},'equality'))
        
            sense_form{1,count} = '=';
            
        elseif(strcmp(constr_type{1,5},'lgr_inequality'))
            
            sense_form{1,count} = '>';
            
        elseif(strcmp(constr_type{1,5},'rgl_inequality'))
            
            sense_form{1,count} = '<';
           
        end % end of if/elseif statements
        
    end
    
    nr_PT2 = (n*m)^2;
    
    for k=1:nr_PT2
        
        count = count + 1;
        
        if(strcmp(constr_type{1,5},'equality'))
        
            sense_form{1,count} = '=';
            
        elseif(strcmp(constr_type{1,5},'lgr_inequality'))
            
            sense_form{1,count} = '>';
            
        elseif(strcmp(constr_type{1,5},'rgl_inequality'))
            
            sense_form{1,count} = '<';
           
        end % end of if/elseif statements
        
    end
    
else
    nr_PT1 = 0;
    nr_PT2 = 0;
end

% McCormick Lifting Inequalities + Box Constraints (Updated 7/10/17)
if(strcmp(set_prop{1,6},'included'))
    
    nr_mcc_lift = 3*(n*m)^(2);
    
    for k=1:nr_mcc_lift
        
        count = count + 1;
        
        if(strcmp(constr_type{1,6},'equality'))
        
            sense_form{1,count} = '=';
            
        elseif(strcmp(constr_type{1,6},'lgr_inequality'))
            
            sense_form{1,count} = '>';
            
        elseif(strcmp(constr_type{1,6},'rgl_inequality'))
            
            sense_form{1,count} = '<';
           
        end % end of if/elseif statements
       
    end
    
else
    
    nr_mcc_lift = 0;
    
end

% Quadratic Lifing Inequalities
if(strcmp(set_prop{1,7},'included'))
    
    % Account for secant inequalities and t_p{ij} upper/lower bounds
    nr_quad_lift = 2*(n*m) + 2*(n*m) + 2*(n*m);

    for k=1:nr_quad_lift
        
        count = count + 1;
        
        if(strcmp(constr_type{1,7},'equality'))
        
            sense_form{1,count} = '=';
            
        elseif(strcmp(constr_type{1,7},'lgr_inequality'))
            
            sense_form{1,count} = '>';
            
        elseif(strcmp(constr_type{1,7},'rgl_inequality'))
            
            sense_form{1,count} = '<';
           
        end % end of if/elseif statements
        
    end
    
else
    nr_quad_lift = 0;
end

if(strcmp(set_prop{1,8},'included'))
    
    nr_box = 2*n*m;
    
    for k=1:nr_box
        
       count = count + 1;
       
       % Either H(i,j) <= Mu(i,j) or -H(i,j) <= -Ld(i,j)
       
       sense_form{1,count} = '<';
       
    end
    
else
    
    nr_box = 0;

end


        
        

