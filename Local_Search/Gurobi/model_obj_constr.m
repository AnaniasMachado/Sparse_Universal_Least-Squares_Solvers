% Objective Function Coefficient vector (model_obj_constr.m) - 6/6/17

function[obj_1n] = model_obj_constr(m,n,set_prop)

% Variables(+order): H (nm), Z (nm), w (2nm), K (nm + nchoosek(nm,2)) 

% Identify if extended formulation (include w/K)
if(strcmp(set_prop{4},'included')|| strcmp(set_prop{5},'included') || strcmp(set_prop{6},'included') || strcmp(set_prop{7},'included'))
    
    % Variables: H (nm), Z (nm), w (2nm), K (nm + nchoosek(nm,2))
    obj_1n = [zeros(1,n*m),ones(1,n*m),zeros(1,2*n*m + n*m + nchoosek(n*m,2))];
    
else
    
    % Variables: H (nm), Z (nm)
    obj_1n = [zeros(1,n*m),ones(1, n*m)];
    
end % end of if/else statement