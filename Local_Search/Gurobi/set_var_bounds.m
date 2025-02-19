% set_var_bounds.m (6/9/17) - In conjunction with nonsym_gurobi.m

% DESCRIPTION: Establish lower bounds on provided variables (or designate
% which are unbounded - in the numerical sense, i.e., less than -1e-20).

% Two possible variable sets:
%
%   (1) LP formulations with H_{ij}, Z_{ij}, K_{ij}{kl}
%       - SVD form of P1, P3, P4, non-symmetric P2, P1 prod. tight, McCormick ineq.
%
%   (2) SOCP formulations with H_{ij}, Z_{ij}, K_{ij}{kl}, w_p{ij}
%       - Quadratic Lifting ineq. (set_prop entry 7)

% INPUT: matrix A dimensions (m,n), constraint set properties (set_prop - char array)

% OUTPUT: var_lb : vector (column) of "unbounded" values -Inf
%         var_ub (optional) : vector (column) of upper bound values


function[var_ub,var_lb] = set_var_bounds(m,n,set_prop)

% Identify which constraint sets are utilized

if(strcmp(set_prop{7},'included'))
    
    % If quadratic lifting inequalities included, then variable space is:
    %   H(i,j), Z(i,j), w_p{ij}, K_{ij}(kl) for all ij (n x m), kl (m x n)
    
    var_lb = -Inf(n*m + n*m + 2*n*m + n*m + nchoosek(n*m,2),1);
    
    % OPTIONAL: Establish a finite upper bound (default in GUROBI is Inf)
    
    var_ub = Inf(n*m + n*m + 2*n*m + n*m + nchoosek(n*m,2),1);
    
elseif((strcmp(set_prop{4},'included') || strcmp(set_prop{5},'included') || strcmp(set_prop{6},'included')) && strcmp(set_prop{7},'not included'))
    
    % Nonsym. P2, P1 prod. tightening, McCormick inequalities utilize the variable space:
    %   H(i,j), Z(i,j), K_{ij}(kl) for all ij (n x m), kl (m x n)
    
    %var_lb = -Inf(n*m + n*m + n*m + nchoosek(n*m,2),1);
    var_lb = -Inf(n*m + n*m + n*m + nchoosek(n*m,2),1);

    
    % OPTIONAL: Establish a finite upper bound
    
    %var_ub = Inf(n*m + n*m + n*m + nchoosek(n*m,2));
    var_ub = Inf(n*m + n*m + n*m + nchoosek(n*m,2),1);
    
else
    
    % SVD of P1, P3, and P4 utilize the variable space: H(i,j), Z(i,j)
    
    var_lb = -Inf(n*m + n*m,1);
    
    % OPTIONAL: Establish an upper bound
    
    var_ub = Inf(n*m + n*m,1);
    
end % end of if/else statement, if(strcmp(set_prop,''))