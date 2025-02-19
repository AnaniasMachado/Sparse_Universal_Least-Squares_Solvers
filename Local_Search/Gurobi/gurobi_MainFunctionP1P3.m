function [norma_gurobi,time_gurobi,H] = gurobi_MainFunctionP1P3(A,m,n)
t1 = tic;

instance = 4;

pre_solve_status = 'enabled';

set_prop = {'included','included','not included','not included','not included','not included','not included','not included','not included','not included','not included'};
%%
[M_s,L,M] = sparse_constr_matrix_gen(A,set_prop,instance);

%%  

[nr,~] = size(M_s);

%%
z_rhs = rhs_constr(m,n,A,set_prop,L,M);
%%

constraint_type = cell(1,8);

% Constraint: P1
constraint_type{1,1} = 'equality';
%constraint_type{1,2} = 'equality';
constraint_type{1,2} = 'lgr_inequality';
constraint_type{1,3} = 'equality';
constraint_type{1,4} = 'equality';
constraint_type{1,5} = 'equality';
constraint_type{1,6} = 'lgr_inequality';
constraint_type{1,7} = 'rgl_inequality';
%%%%%constraint_type{1,8} = 'equality';
constraint_type{1,8} = 'lgr_inequality';
constraint_type{1,9} = 'equality';

objective_type = cell(1,1);

objective_type{1,1} = 'min';

%%

[sense_form] = constr_sense(m,n,A,nr,constraint_type,set_prop);

%%

obj_1n = model_obj_constr(m,n,set_prop);

%%

if(strcmp(pre_solve_status,'disabled'))
    
    params.predeprow = -1;
    
    params.presolve = 0;
    
elseif(strcmp(pre_solve_status,'enabled'))
    
    params.predeprow = -1;
    
    params.presolve = -1;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GUROBI MODEL CONSTRUCTION:

[var_ub,var_lb] = set_var_bounds(m,n,set_prop);

model.A = sparse(M_s);
model.obj = obj_1n;
model.modelsense = 'Min';
model.rhs = z_rhs;
%keyboard
model.sense = char(sense_form);
model.lb = var_lb;
%%% model.lb = -Inf*ones(m*n*2,1); Esse era o problema de antes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result = gurobi(model,params);

% Check to see if problem is infeasible or not; then display obj. value
if(isfield(result, 'objval'))
    
    norma_gurobi = result.objval;
    
end % end of if statment, if(isfield(result, 'objval'))

% If problem feasible (and opt. soln found), extract variable information
if(isfield(result, 'x'))

    [H,Z,MKIJ] = extract_variables(m,n,result.x,set_prop);

end
%keyboard
H = sparse(H);


time_gurobi = toc(t1);
end
