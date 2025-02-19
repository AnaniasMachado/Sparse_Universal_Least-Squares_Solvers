
function[M_sparse,L,M] = sparse_constr_matrix_gen(A,set_prop,instance)

% Number of rows/columns of input A
[m,n] = size(A);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ||H||_{1} : -Z <= H <= Z equivalent formulation Z - H >= 0 and H + Z >= 0

% H + Z >= 0
one_NL_row_index = [];
one_NL_col_index = [];
one_NL_mat_val = [];
one_NL_sparse = [];

% Z - H >= 0
one_NU_row_index = [];
one_NU_col_index = [];
one_NU_mat_val = [];
one_NU_sparse = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P1_row_index = [];
P1_col_index = [];
P1_mat_val = [];
P1_sparse = [];

P3_row_index = [];
P3_col_index = [];
P3_mat_val = [];
P3_sparse = [];

P4_row_index = [];
P4_col_index = [];
P4_mat_val = [];
P4_sparse = [];

%%% Property 10: Symmetric

P_sym_row_index = [];
P_sym_col_index = [];
P_sym_mat_val = [];
P_sym_sparse = [];

%%% Property 11: P2 Linearized (To be used together with P1 and P3)

P2L_row_index = [];
P2L_col_index = [];
P2L_mat_val = [];
P2L_sparse = [];

%%%%%%%%%%%%%%%%%%%%%%

PT1_row_index = [];
PT1_col_index = [];
PT1_mat_val = [];
PT1_sparse = [];

PT2_row_index = [];
PT2_col_index = [];
PT2_mat_val = [];
PT2_sparse = [];

%%%%%%%%%%%%%%%%%%%%%%

NS_row_index = [];
NS_col_index = [];
NS_mat_val = [];
NS_sparse = [];

mccl_row_index = [];
mccl_col_index = [];
mccl_mat_val = [];
mccl_sparse = [];

% NOTE(7/11/17): Box constraints for McCormick

% Lower bounds : Lambda (Ld), Ld <= H
Ld_row_index = [];
Ld_col_index = [];
Ld_mat_val = [];
Ld_sparse = [];

% Upper bounds : Mu (Mu), H <= Mu
Mu_row_index = [];
Mu_col_index = [];
Mu_mat_val = [];
Mu_sparse = [];

%%%%%%%%%%%%%%%%%%%%%%%

% Secant constraints (h_{ij},w_p{ij} variables)
secant_row_index = [];
secant_col_index = [];
secant_mat_val = [];
secant_sparse = [];

% Bounds on h_{ij} (upper and lower bounds)
tb_row_index = [];
tb_col_index = [];
tb_mat_val = [];
tb_sparse = [];

ta_row_index = [];
ta_col_index = [];
ta_mat_val = [];
ta_sparse = [];

% Quadratic lifting ineqs (h_{ij}, \mathcal{K}_{ij}, w_p{ij} variables)
Qc_row_index = [];
Qc_col_index = [];
Qc_mat_val = [];
Qc_sparse = [];

q_row_index = [];
q_col_index = [];
q_mat_val = [];
q_sparse = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Box constraints on H : L(i,j) <= H(i,j) <= M(i,j)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "data_type" denotes if bounds generated of known data/matrix
%   (type 1) : "provided" -> use data_mat
%   (type 2) : "generate" -> data_mat = [], generate via another procedure

data_type = 'provided';

data_source = 'mpp'; % Added 7/26/17

H_box = pinv(A);

[L,M] = box_constr_gen(m,n,'provided',H_box);


for i=1:n
    for j=1:m
        
        ij = (i-1)*m + j;
        
        % Constraint: H(i,j) + Z(i,j) >=0
        
        % H(i,j) coefficient
        one_NL_row_index = [one_NL_row_index, ij];
        one_NL_col_index = [one_NL_col_index, ij]; % H(i,j) before Z(i,j)
        one_NL_mat_val = [one_NL_mat_val, 1];
        
        % Z(i,j) coefficient
        one_NL_row_index = [one_NL_row_index, ij];
        one_NL_col_index = [one_NL_col_index, (n*m) + ij];
        one_NL_mat_val = [one_NL_mat_val, 1];
        
    end
end


one_NL_sparse = sparse(one_NL_row_index, one_NL_col_index, one_NL_mat_val, n*m, n*m + n*m);;


% Inequalities: Z - H >= 0 (Rewritten from H <= Z)

for i=1:n
    for j=1:m
        
        ij = (i-1)*m + j;
        
        % Constraint: Z(i,j) - H(i,j) >= 0
        
        % Z(i,j) coefficient
        one_NU_row_index = [one_NU_row_index, ij];
        one_NU_col_index = [one_NU_col_index, (n*m) + ij];
        one_NU_mat_val = [one_NU_mat_val, 1];
        
        % H(i,j) coefficient
        one_NU_row_index = [one_NU_row_index, ij];
        one_NU_col_index = [one_NU_col_index, ij];
        one_NU_mat_val = [one_NU_mat_val, -1];
        
    end
end

one_NU_sparse = sparse(one_NU_row_index, one_NU_col_index, one_NU_mat_val, n*m, n*m + n*m);

% Identify which MPP properties are to be included

% P1: AHA == A

% SVD form: (S*V')*H*(U*S) == U'*A*V,
%   Dimesions: U (m x r), S (r x r), V (n x r)

[U,S,V] = svd(A); % NOTE: A = U*S*V'

r = rank(A);

% U : reduce number of columns (m x m) -> (m x r)

U = U(:,1:r);

% S : reduce number of columns and rows (m x n) -> (r x r)

S = S(1:r,1:r);

% V : reduce number of columns (n x n) -> (n x r)

V = V(:,1:r);

for i=1:r
    for j=1:r
        
        ij = (i-1)*r + j;
        
        for l=1:n
            for k=1:m
                
                lk = (l-1)*m + k;
                
                P1_row_index = [P1_row_index,ij];
                P1_col_index = [P1_col_index,lk];
                
                % h_{lk} coefficient
                coeff_hlk = (S(i,i)^(2))*V(l,i)*U(k,j);
                
                P1_mat_val = [P1_mat_val,coeff_hlk];
                
            end % end of for loop, for k=1:M
        end % end of for loop, for l=1:n
        
    end % end of for loop, for j=1:n
end % end of for loop, for i=1:m



% Check for P1-tightening, Non-sym P2, McCormick Ineq, or Quadratic Ineq


% Exclude w_p{ij}, \mathcal{K} variables, simplify LP formulation
P1_sparse = sparse(P1_row_index,P1_col_index,P1_mat_val,r*r,n*m + n*m);

% Property 3: (AH)' = AH <=> (AH)' - AH = 0
if(strcmp(set_prop{2},'included'))
    
    count = 0; % Row index
    
    for i=1:m
        for j=1:m
            
            %ij = (i-1)*n + j; % mapping index - by row of (AH)/(AH)'
            
            if(i ~= j && j > i)
                
                count = count + 1;
                
                for k=1:n
                    
                    % Variable h_{kj} : (AH)
                    
                    kj = (k-1)*m + j;
                    
                    P3_row_index = [P3_row_index,count];
                    P3_col_index = [P3_col_index,kj];
                    P3_mat_val = [P3_mat_val,-A(i,k)];
                    
                    % Variable h_{ki} : (AH)'
                    
                    ki = (k-1)*m + i;
                    
                    P3_row_index = [P3_row_index,count];
                    P3_col_index = [P3_col_index,ki];
                    P3_mat_val = [P3_mat_val,A(j,k)];
                    
                end % end of for loop, for k=1:n
                
            end % end of if statement, if(i ~= j)
            
        end % end of inner for loop, for j=1:n
    end % end of outer for loop, for i=1:m
    
    % Check for P1-tightening, Non-sym P2, McCormick Ineq, or Quadratic Ineq
    if(strcmp(set_prop{4},'included') || strcmp(set_prop{5},'included') || strcmp(set_prop{6},'included') || strcmp(set_prop{7},'included'))
        
        % UPDATE (6/6/17): Need to account for Z (nm), w (2nm) variables
        
        P3_sparse = sparse(P3_row_index,P3_col_index,P3_mat_val,count,(n*m + n*m + 2*n*m) + n*m + nchoosek(n*m,2));
        
    else
        
        P3_sparse = sparse(P3_row_index,P3_col_index,P3_mat_val,count,n*m + n*m);
        
    end % end of if/else statement
    
else
    
    P3_sparse = [];
    
end % end of if/else statement, if(strcmp(set_prop{2},''))


% Property 4: (HA)' = HA <=> (HA)' - HA = 0
if(strcmp(set_prop{3},'included'))
    
    count = 0; % Row index
    
    for i=1:n
        for j=1:n
            
            %ij = (i-1)*m + j; % mapping index - by row of (HA)/(HA)'
            
            if(i ~= j && j > i)
                
                count = count + 1;
                
                for k=1:m
                    
                    % Variable h_{ik} : (HA)
                    
                    ik = (i-1)*m + k;
                    
                    P4_row_index = [P4_row_index,count];
                    P4_col_index = [P4_col_index,ik];
                    P4_mat_val = [P4_mat_val,-A(k,j)];
                    
                    % Variable h_{jk} : (HA)'
                    
                    jk = (j-1)*m + k;
                    
                    P4_row_index = [P4_row_index,count];
                    P4_col_index = [P4_col_index,jk];
                    P4_mat_val = [P4_mat_val,A(k,i)];
                    
                end % end of for loop, for k=1:n
                
            end % end of if statement, if(i ~= j)
            
        end % end of inner for loop, for j=1:n
    end % end of outer for loop, for i=1:m
    
    % Check for P1-tightening, Non-sym P2, McCormick Ineq, or Quadratic Ineq
    if(strcmp(set_prop{4},'included') || strcmp(set_prop{5},'included') || strcmp(set_prop{6},'included') || strcmp(set_prop{7},'included'))
        
        % UPDATE(6/6/17): Incorporate Z (nm) and w (2nm) variables
        
        P4_sparse = sparse(P4_row_index,P4_col_index,P4_mat_val,count,(n*m + n*m + 2*n*m) + n*m + nchoosek(n*m,2));
        
    else
        
        P4_sparse = sparse(P4_row_index,P4_col_index,P4_mat_val,count,n*m + n*m);
        
    end % end of if/else statement
    
else
    
    P4_sparse = [];
    
end % end of if/else statement, if(strcmp(set_prop{3},''))

% Property 10: Equalities: H' - H == 0 , for m=n
if strcmp(set_prop{10},'included')
    count_sym=0;
    for i=1:m
        for j= i+1 :m
            count_sym = count_sym + 1;
            ij = (i-1)*m + j;
            ji=(j-1)*m+i;
            % Constraint: H(j,i) - H(i,j) = 0
            
            % H(i,j) coefficient
            P_sym_row_index = [P_sym_row_index, count_sym];
            P_sym_col_index = [P_sym_col_index, ij];
            P_sym_mat_val = [P_sym_mat_val, 1];
            
            % H(j,i) coefficient
            P_sym_row_index = [P_sym_row_index, count_sym];
            P_sym_col_index = [P_sym_col_index, ji];
            P_sym_mat_val = [P_sym_mat_val, -1];
            
            
            count_sym = count_sym + 1;
            ij = (i-1)*m + j;
            ji=(j-1)*m+i;
            % Constraint: H(j,i) - H(i,j) = 0
            
            % H(i,j) coefficient
            P_sym_row_index = [P_sym_row_index, count_sym];
            P_sym_col_index = [P_sym_col_index, ij];
            P_sym_mat_val = [P_sym_mat_val, -1];
            
            % H(j,i) coefficient
            P_sym_row_index = [P_sym_row_index, count_sym];
            P_sym_col_index = [P_sym_col_index, ji];
            P_sym_mat_val = [P_sym_mat_val, 1];
        end
    end
    
%    keyboard
    P_sym_sparse = sparse(P_sym_row_index, P_sym_col_index, P_sym_mat_val, count_sym, 2*n*m);
    
else
    
    P_sym_sparse = [];
    
end


% Property 2 Linearized: HAA+ = H <=> HAA+ - H = 0
if(strcmp(set_prop{11},'included'))
    Mat = A*pinv(A);
    count = 0; % Row index
    
    for i=1:n
        for j=1:m
            
            ij = (i-1)*m + j; 
            count = count + 1;
            
            P2L_row_index = [P2L_row_index,count];
            P2L_col_index = [P2L_col_index,ij];
            P2L_mat_val = [P2L_mat_val,-1];
                
            for k=1:m                
                    % Variable h_{ik} : (HMat)
                    
                    ik = (i-1)*m + k;
                    
                    P2L_row_index = [P2L_row_index,count];
                    P2L_col_index = [P2L_col_index,ik];
                    P2L_mat_val = [P2L_mat_val,Mat(k,j)];
             end % end of for loop, for k=1:m
            
        end % end of inner for loop, for j=1:m
    end % end of outer for loop, for i=1:n
        
%    keyboard
    P2L_sparse = sparse(P2L_row_index, P2L_col_index, P2L_mat_val, count, 2*n*m);
    
else
    
    P2L_sparse = [];
    
end

% Nonsymmetric P2 (Linear in \mathcal{K}_{ji}): HAH - H = 0
if(strcmp(set_prop{4},'included'))
    
    for j=1:n
        for i=1:m
            
            ji = (j-1)*m + i; % mapping index
            
            % Sparse format
            
            NS_row_index = [NS_row_index,ji];
            NS_col_index = [NS_col_index,ji];
            
            NS_mat_val = [NS_mat_val,-1];
            
            % Dense format
            %NS(ji,ji) = -1;
            
            for k=1:n
                for l=1:m
                    
                    % K[ji](lk) := h_{jl}h_{ki}
                    
                    lk = (k-1)*m + l;
                    
                    if(ji <= lk)
                        
                        t_sum = sum(1:ji-1);
                        
                        NS_row_index = [NS_row_index,ji];
                        
                        % UPDATE(6/6/17): Incorporate Z (nm) and w (2nm)
                        %   variables in index (order: H,Z,w,K) - before just H,K
                        NS_col_index = [NS_col_index,n*m + (2*n*m) + (ji)*(n*m) + lk - t_sum];
                        
                        NS_mat_val = [NS_mat_val,A(l,k)];
                        
                        %NS(ji, (jl)*(n*m) + ki - t_sum) = A(l,k);
                        
                    elseif(ji > lk)
                        
                        % UPDATE(8/3/17): Correcting error in generation
                        
                        % Use K[kl](ij) instead of K[ji](lk)
                        
                        kl = (k-1)*m + l; % Index for which matrix (nm total)
                        
                        ij = (j-1)*m + i; % Index for entry in a (m x n) matrix
                        
                        t_sum = sum(1:kl-1);
                        
                        NS_row_index = [NS_row_index,ji];
                        
                        NS_col_index = [NS_col_index,n*m + (2*n*m) + (kl)*(n*m) + ij - t_sum];
                        
                        NS_mat_val = [NS_mat_val,A(l,k)];
                        
                    end % end of if statement, if(jl <= ki)
                    
                end % end of inner for loop, for l=1:m
            end % end of outer for loop, for k=1:n
            
        end % end of inner for loop, for i=1:m
    end % end of outer for loop, for j=1:n
    
    % UPDATE(6/6/17): Incorporate Z (nm), w (2nm) variables
    
    NS_sparse = sparse(NS_row_index,NS_col_index,NS_mat_val,n*m,(n*m + n*m + 2*n*m) + n*m + nchoosek(n*m,2));
    
else
    
    NS_sparse = [];
    
end % end of if statement, if(strcmp(set_prop(4),'')) - nonsym P2


% Product Tightening of P1 (Linear in \mathcal{K}_{ji} and h_{ji})

if(strcmp(set_prop{5},'included'))
    
    % (i) <A,K[j,i]> = H(j,i) equivalent to sum(k=1:m)sum(l=1:n)K[j,i](k,l)A(k,l)
    
    for j=1:n
        for i=1:m
            
            ji = (j-1)*m + i; % Mapping index for K[j,i] (which one)
            
            % Constraint coefficient for entry H(j,i):
            PT1_row_index = [PT1_row_index,ji];
            
            PT1_col_index = [PT1_col_index,ji];
            
            PT1_mat_val = [PT1_mat_val,-1];
            
            % For each entry of H, H(j,i): construct left-hand sum...
            %
            %   sum(k=1:m) sum(l=1:n) K[ji](k,l)A(k,l) for all ji pairs
            
            for k=1:m
                for l=1:n
                    
                    kl = (l-1)*m + k; % Mapping index for ENTRY (k,l) in K[j,i]
                    
                    % Sum of all entries in K[j,i], K[ji](k,l), multiplied by
                    % each entry in input matrix A, A(k,l), for all K[j,i]
                    
                    PT1_row_index = [PT1_row_index,ji];
                    
                    if(ji <= kl)
                        
                        t_sum = sum(1:ji-1);
                        
                        % Index for K[j,i](k,l)
                        PT1_col_index = [PT1_col_index,n*m + (2*n*m) + (ji)*(n*m) + kl - t_sum];
                        
                    elseif(kl < ji)
                        
                        t_sum = sum(1:kl-1);
                        
                        % Index for K[l,k](i,j) [equivalent to K[j,i](k,l)]
                        PT1_col_index = [PT1_col_index,n*m + (2*n*m) + (kl)*(n*m) + ji - t_sum];
                        
                    end % end of if/elseif statement, if(ji <= kl || kl > ji)
                    
                    % Corresponding coefficient A(k,l) for K[j,i](k,l)
                    PT1_mat_val = [PT1_mat_val,A(k,l)];
                    
                    
                end % end of for loop, for l=1:n
            end % end of for loop, for k=1:m
            
        end % end of for loop, for i=1:m
    end % end of for loop, for j=1:n
    
    % Notation: \mathcal{K}_{pk}{ql} := H(p,q)H(l,k)
    % - p=1:n, q=1:m, k=1:m, l=1:n
    
    % (ii) sum(l=1:n)sum(k=1:m)A(i,l)A(k,j)K[i,j](k,l) = A(i,j)H(p,q) for all ij (n x m), pq (m x n)
    
    for p=1:n
        for k=1:m
            for j=1:n
                for i=1:m
                    
                    ji = (j-1)*m + i; % Mapping index (which A(i,j))
                    
                    pk = (p-1)*m + k; % Mapping index (which \mathcal{K}_{pk})
                    
                    for q=1:m
                        for l=1:n
                            
                            ql = (l-1)*m + q; % Mapping index (which entry)
                            
                            pq = (p-1)*m + q; % Mapping index (establish row index)
                            
                            if(pk <= ql)
                                
                                % Sparse format :  -A(i,j)
                                PT2_row_index = [PT2_row_index, ji + (pq-1)*(n*m)];
                                
                                PT2_col_index = [PT2_col_index,pq];
                                
                                PT2_mat_val = [PT2_mat_val,-A(i,j)];
                                
                                t_sum = sum(1:pk-1);
                                
                                % Sparse format : A(i,l)*A(k,j)
                                PT2_row_index = [PT2_row_index,ji + (pq-1)*(n*m)];
                                
                                % Variable order: H,Z,w,K (UPDATE[6/6/17])
                                PT2_col_index = [PT2_col_index,n*m + (2*n*m) + (pk)*(n*m) + ql - t_sum];
                                
                                PT2_mat_val = [PT2_mat_val,A(i,l)*A(k,j)];
                                
                            elseif(ql < pk)
                                
                                % Sparse format :  -A(i,j)
                                PT2_row_index = [PT2_row_index, ji + (pq-1)*(n*m)];
                                
                                PT2_col_index = [PT2_col_index,pq];
                                
                                PT2_mat_val = [PT2_mat_val,-A(i,j)];
                                
                                t_sum = sum(1:ql-1);
                                
                                % Sparse format : A(i,l)*A(k,j)
                                PT2_row_index = [PT2_row_index,ji + (pq-1)*(n*m)];
                                
                                % UPDATE(6/6/17): variables H,Z,w,K
                                PT2_col_index = [PT2_col_index,n*m + 2*n*m + (ql)*(n*m) + pk - t_sum];
                                
                                PT2_mat_val = [PT2_mat_val,A(i,l)*A(k,j)];
                                
                            end % end of if statement
                        end % end of inner for loop, for q=1:m
                    end % end of outer for loop, for l=1:n
                    
                end % end of inner for loop, for i=1:m
            end % end of outer for loop, for j=1:n
            
        end % end of inner for loop, for q=1:m
    end % end of outer for loop, for p=1:n
    
    % UPDATE(6/6/17): Incorporate Z (nm), w (2nm) variables
    
    PT1_sparse = sparse(PT1_row_index,PT1_col_index,PT1_mat_val,n*m,(n*m + n*m + 2*n*m) + n*m + nchoosek(n*m,2));
    
    PT2_sparse = sparse(PT2_row_index,PT2_col_index,PT2_mat_val,(n*m)^2,(n*m + n*m + 2*n*m) + n*m + nchoosek(n*m,2));
    
    PT_sparse = [PT1_sparse; PT2_sparse];
    
else
    
    PT_sparse = [];
    
end % end if/else statement, if(strcmp(set_prop(5),''))

% McCormick Lifting Inequalities
if(strcmp(set_prop{6},'included'))
    
    % Box constraints on H : L(i,j) <= H(i,j) <= M(i,j)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % "data_type" denotes if bounds generated of known data/matrix
    %   (type 1) : "provided" -> use data_mat
    %   (type 2) : "generate" -> data_mat = [], generate via another procedure
    
    data_type = 'provided';
    
    data_source = 'mpp'; % Added 7/26/17
    
    if(strcmp(data_type,'provided'))
        
        if(strcmp(data_source,'existing'))
            
            data_load_spec = 'SPIg_%d.mat';
            
            data_load = sprintf(data_load_spec,instance);
            
            load(data_load,'H');
            
            H_box = H;
            
        elseif(strcmp(data_source,'mpp'))
            
            H_box = pinv(A);
            
        end % end of if/elseif statement
        
        [L,M] = box_constr_gen(m,n,'provided',H_box);
        
    end % end of if statement
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    count = 0;
    
    for i=1:n
        for j=1:m
            
            ij = (i-1)*m + j;
            
            for k=1:m
                for l=1:n
                    
                    % For h_{ik}, h_{lj} : mapping indices
                    
                    ik = (i-1)*m + k; % h_{ik}
                    lj = (l-1)*m + j; % h_{lj}
                    
                    kl = (l-1)*m + k;
                    
                    % McCormick Inequalities : Type 1-3 for h_{ik}, h_{lj}
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % TYPE 1 :
                    
                    count = count + 1;
                    
                    % h_{ik}
                    
                    mccl_row_index = [mccl_row_index,count];
                    mccl_col_index = [mccl_col_index,ik];
                    mccl_mat_val = [mccl_mat_val,-L(l,j)];
                    
                    % h_{lj}
                    
                    mccl_row_index = [mccl_row_index,count];
                    mccl_col_index = [mccl_col_index,lj];
                    mccl_mat_val = [mccl_mat_val,-L(i,k)];
                    
                    
                    % \mathcal{K}_{ij}{kl}
                    
                    mccl_row_index = [mccl_row_index,count];
                    
                    if(ij <= kl)
                        
                        sum_t = sum(1:ij-1);
                        
                        % UPDATE(6/6/17): Incorporate Z (nm), w (2nm)
                        
                        mccl_col_index = [mccl_col_index,n*m + (2*n*m) + (ij)*n*m + kl - sum_t];
                        mccl_mat_val = [mccl_mat_val,1];
                        
                    elseif(kl < ij)
                        
                        sum_t = sum(1:kl-1);
                        
                        % UPDATE(6/6/17): Incorporate Z (nm), w (2nm)
                        
                        mccl_col_index = [mccl_col_index,n*m + 2*n*m + (kl)*n*m + ij - sum_t];
                        mccl_mat_val = [mccl_mat_val,1];
                        
                    end % end of if/elseif statement
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % TYPE 2 :
                    
                    count = count + 1;
                    
                    % h_{ik}
                    
                    mccl_row_index = [mccl_row_index,count];
                    mccl_col_index = [mccl_col_index,ik];
                    mccl_mat_val = [mccl_mat_val,M(l,j)];
                    
                    % h_{lj}
                    
                    mccl_row_index = [mccl_row_index,count];
                    mccl_col_index = [mccl_col_index,lj];
                    mccl_mat_val = [mccl_mat_val,M(i,k)];
                    
                    
                    % \mathcal{K}_{ij}{kl}
                    
                    mccl_row_index = [mccl_row_index,count];
                    
                    if(ij <= kl)
                        
                        sum_t = sum(1:ij-1);
                        
                        % UPDATE(6/6/17): Incorporate Z (nm), w (2nm)
                        
                        mccl_col_index = [mccl_col_index,n*m + 2*n*m + (ij)*n*m + kl - sum_t];
                        mccl_mat_val = [mccl_mat_val,-1];
                        
                    elseif(kl < ij)
                        
                        sum_t = sum(1:kl-1);
                        
                        % UPDATE(6/6/17): Incorporate Z (nm), w (2nm)
                        
                        mccl_col_index = [mccl_col_index,n*m + 2*n*m + (kl)*n*m + ij - sum_t];
                        mccl_mat_val = [mccl_mat_val,-1];
                        
                    end % end of if/elseif statement
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % TYPE 3 :
                    
                    count = count + 1;
                    
                    % h_{ik}
                    
                    mccl_row_index = [mccl_row_index,count];
                    mccl_col_index = [mccl_col_index,ik];
                    mccl_mat_val = [mccl_mat_val,-M(l,j)];
                    
                    % h_{lj}
                    
                    mccl_row_index = [mccl_row_index,count];
                    mccl_col_index = [mccl_col_index,lj];
                    mccl_mat_val = [mccl_mat_val,-M(i,k)];
                    
                    
                    % \mathcal{K}_{ij}{kl}
                    
                    mccl_row_index = [mccl_row_index,count];
                    
                    if(ij <= kl)
                        
                        sum_t = sum(1:ij-1);
                        
                        % UPDATE(6/6/17): Incorporate Z (nm), w (2nm)
                        
                        mccl_col_index = [mccl_col_index,n*m + 2*n*m + (ij)*n*m + kl - sum_t];
                        mccl_mat_val = [mccl_mat_val,1];
                        
                    elseif(kl < ij)
                        
                        sum_t = sum(1:kl-1);
                        
                        % UPDATE(6/6/17): Incorporate Z (nm), w (2nm)
                        
                        mccl_col_index = [mccl_col_index,n*m + 2*n*m + (kl)*n*m + ij - sum_t];
                        mccl_mat_val = [mccl_mat_val,1];
                        
                    end % end of if/elseif statement
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                end % end of for loop, for l=1:n
            end % end of for loop, for k=1:m
            
        end % end of for loop, for j=1:m
    end % end of for loop, for i=1:n
    
    mccl_sparse = sparse(mccl_row_index,mccl_col_index,mccl_mat_val,3*(n*m)^2,(n*m + n*m + 2*n*m) + n*m + nchoosek(n*m,2));
    
else
    
    mccl_sparse = [];
    
end % end of if statement, if(strcmp(set_prop(6),''))

% Quadratic Lifting Inequalities
if(strcmp(set_prop{7},'included'))
    
    % Provided data (L,M,U,V)
    %   [1] L(i,j) <= H <= M(i,j) for all (ij) in (nxm)
    %   [2] U^(ij) mx1 vectors, V^(ij) nx1 vectors for all (ij) in (nxm)
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Generate a_k{ij} (alpha) and b_k{ij} (beta) for k=1,2 and (ij) in (nxm)
    
    a = zeros(2,n*m); % a(1,:) a_1{ij}, a(2,:) a_2{ij} for all ij
    b = zerso(2,n*m); % b(1,:) b_1{ij}, b(2,:) b_2{ij} for all ij
    
    for p=1:2
        
        % Generate all a_k{ij}/b_k{ij}
        
        for i=1:n
            for j=1:m
                
                ij = (i-1)*m + j;
                
                sum_am = 0;
                sum_an = 0;
                sum_bm = 0;
                sum_bn = 0;
                
                for l=1:m
                    
                    % Holds for both k=1 and k=2
                    
                    min_am = min(U(l,ij)*L(i,l),U(l,ij)*M(i,l));
                    max_bm = max(U(l,ij)*L(i,l),U(l,ij)*M(i,l));
                    
                    sum_am = sum_am + min_am;
                    sum_bm = sum_bm + max_bm;
                    
                end % end of for loop, for l=1:m
                
                for l=1:n
                    
                    min_an = min(V(l,ij)*L(l,j),V(l,ij)*M(l,j));
                    max_bn = max(V(l,ij)*L(l,j),V(l,ij)*M(l,j));
                    
                    sum_an = sum_an + min_an;
                    sum_bn = sum_bn + max_bn;
                    
                end % end of for loop, for l=1:n
                
                
                a(p,ij) = (1/2)*(sum_am + ((-1)^(p-1))*sum_an);
                b(p,ij) = (1/2)*(sum_bm + ((-1)^(p-1))*sum_bn);
                
            end % end of for loop, for j=1:m
        end % end of for loop, for i=1:n
        
    end % end of for loop, for k=1:2
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Notation: H matrix variable, U/V given
    
    % t_1{ij} := (1/2)(U(:,ij)'H(i,:)' + V(:,ij)'H(:,j))
    % t_2{ij} := (1/2)(U(:,ij)'H(i,:)' - V(:,ij)'H(:,j))
    
    % Secant form: w_p{ij} scalar variables, a/b depend on U/V given
    
    % -[(t_p{ij} - a(p,ij))(b(p,ij)^{2} - a(p,ij)^{2})/(b(p,ij) - a(p,ij)) + a(p,ij)^{2} <= w(p,ij)
    
    % Secant constraints:
    
    count_sec = 0; % Track # of rows generated
    
    for p=1:2
        
        for i=1:n
            for j=1:m
                
                ij = (i-1)*m + j; % mapping index, a(p,ij)/b(p,ij)
                
                pij = (p-1)*n*m + ij; % Row index for secant constraints
                
                count_sec = count_sec + 1;
                
                % Coefficients for h_{ik} (k=1:m) and h_{kj} (k=1:n)
                for k=1:max(m,n)
                    
                    if(k <= m && max(m,n) == n)
                        
                        % For h_{ik}
                        ik = (i-1)*m + k;
                        
                        % For h_{kj}
                        kj = (k-1)*m + j;
                        
                        if(ik == kj)
                            
                            secant_row_index = [secant_row_index,count_sec];
                            secant_col_index = [secant_col_index,ik]; % or kj
                            
                            coeff_val = -(1/2)*((b(p,ij)^2 - a(p,ij)^2)/(b(p,ij) - a(p,ij)))*(U(k,ij) + ((-1)^(p-1))*V(k,ij));
                            
                            secant_mat_val = [secant_mat_val,coeff_val]
                            
                        elseif(ik ~= kj)
                            
                            % h_{ik} coefficient term
                            secant_row_index = [secant_row_index,count_sec];
                            secant_col_index = [secant_col_index,ik];
                            
                            coeff_hik = -(1/2)*((b(p,ij)^2 - a(p,ij)^2)/(b(p,ij) - a(p,ij)))*U(k,ij);
                            
                            secant_mat_val = [secant_mat_val,coeff_hik];
                            
                            % h_{kj} coefficient term
                            secant_row_index = [secant_row_index,count_sec];
                            secant_col_index = [secant_col_index,kj];
                            
                            coeff_hkj = -(1/2)*((b(p,ij)^2 - a(p,ij)^2)/(b(p,ij) - a(p,ij)))*((-1)^(p-1));
                            
                            secant_mat_val = [secant_mat_val,coeff_hkj];
                            
                        end % end of if/elseif statement, if(ik == kj)
                        
                    else
                        
                        kj = (k-1)*m + j; % k=1:n
                        
                        % k > m, only generate h_{kj}
                        secant_row_index = [secant_row_index,count_sec];
                        secant_col_index = [secant_col_index,kj];
                        
                        coeff_hkj = -(1/2)*((b(p,ij)^2 - a(p,ij)^2)/(b(p,ij) - a(p,ij)))*((-1)^(p-1));
                        
                        secant_mat_val = [secant_mat_val,coeff_hkj];
                        
                    end % end of if/else statement, if(k <= m && max(m,n) == n)
                    
                    if(k <= n && max(m,n) == m)
                        
                        % For h_{ik}
                        ik = (i-1)*m + k;
                        
                        % For h_{kj}
                        kj = (k-1)*m + j;
                        
                        if(ik == kj)
                            
                            secant_row_index = [secant_row_index,count_sec];
                            secant_col_index = [secant_col_index,ik]; % or kj
                            
                            coeff_val = -(1/2)*((b(p,ij)^2 - a(p,ij)^2)/(b(p,ij) - a(p,ij)))*(U(k,ij) + ((-1)^(p-1))*V(k,ij));
                            
                            secant_mat_val = [secant_mat_val,coeff_val]
                            
                        elseif(ik ~= kj)
                            
                            % h_{ik} coefficient term
                            secant_row_index = [secant_row_index,count_sec];
                            secant_col_index = [secant_col_index,ik];
                            
                            coeff_hik = -(1/2)*((b(p,ij)^2 - a(p,ij)^2)/(b(p,ij) - a(p,ij)))*U(k,ij);
                            
                            secant_mat_val = [secant_mat_val,coeff_hik];
                            
                            % h_{kj} coefficient term
                            secant_row_index = [secant_row_index,count_sec];
                            secant_col_index = [secant_col_index,kj];
                            
                            coeff_hkj = -(1/2)*((b(p,ij)^2 - a(p,ij)^2)/(b(p,ij) - a(p,ij)))*((-1)^(p-1));
                            
                            secant_mat_val = [secant_mat_val,coeff_hkj];
                            
                        end % end of if/elseif statement, if(ik == kj)
                        
                    else
                        
                        ik = (i-1)*m + k; % k=1:m
                        
                        % k > n, only generate h_{ik}
                        secant_row_index = [secant_row_index,count_sec];
                        secant_col_index = [secant_col_index,ik];
                        
                        coeff_hik = -(1/2)*((b(p,ij)^2 - a(p,ij)^2)/(b(p,ij) - a(p,ij)))*((-1)^(p-1));
                        
                        secant_mat_val = [secant_mat_val,coeff_hik];
                        
                    end % end of if/else statement, if(k <= n && max(m,n) == m)
                    
                end % end of for loop, for k=1:max(m,n)
                
                % Coefficient for w_p{ij} variables
                
                secant_row_index = [secant_row_index,count_sec];
                
                % UPDATE(6/6/17): Account for H(i,j) and Z(i,j), shift by 2*n*m
                secant_col_index = [secant_col_index,n*m + n*m + (p-1)*nm + ij];
                secant_mat_val = [secant_mat_val,-1]; % Move w_p{ij} term to left side of inequality
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Right-hand side (reference)
                
                % Form: a(p,ij)^2 - ((b(p,ij)^2 - a(p,ij)^2)/(b(p,ij) - a(p,ij)))*a(p,ij)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            end % end of for loop, for j=1:m
        end % end of for loop, for i=1:n
        
    end % end of for loop, for p=1:2
    
    % Inequalities t_p{ij} <= b_p{ij} and -t_p{ij} <= -a_p{ij}
    
    % t_1{ij} := (1/2)(U(:,ij)'H(i,:)' + V(:,ij)'H(:,j))
    % t_2{ij} := (1/2)(U(:,ij)'H(i,:)' - V(:,ij)'H(:,j))
    
    counta = 0;
    countb = 0;
    
    for p=1:2
        
        for i=1:n
            for j=1:m
                
                ij = (i-1)*m + j;
                
                pij = (p-1)*n*m + ij; % An individual # constraint count
                
                % Track # of constraints generated
                counta = counta + 1;
                countb = countb + 1;
                
                for k=1:max(m,n)
                    
                    if(k <= m && max(m,n) == n)
                        
                        % For h_{ik}
                        ik = (i-1)*m + k; % k=1:m
                        
                        % For h_{kj}
                        kj = (k-1)*m + j; % k=1:n
                        
                        if(ik == kj)
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % Inequality: t_p{ij} <= b(p,ij)
                            
                            % Coefficients for h_{ik} and h_{kj} combined
                            tb_row_index = [tb_row_index,countb];
                            tb_col_index = [tb_col_index,ik]; % or kj
                            
                            coeff_val_b = (1/2)*(U(k,ij) + ((-1)^(p-1))*V(k,ij));
                            
                            tb_mat_val = [tb_mat_val,coeff_val_b];
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % Inequality: -t_p{ij} <= -a(p,ij)
                            
                            % Coefficients for h_{ik} and h_{kj} combined
                            ta_row_index = [ta_row_index,counta];
                            ta_col_index = [ta_col_index,ik]; % or kj
                            
                            coeff_val_a = -(1/2)*(U(k,ij) + ((-1)^(p-1))*V(k,ij));
                            
                            ta_mat_val = [ta_mat_val,coeff_val_a];
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                            
                        elseif(ik ~= kj)
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % Inequality: t_p{ij} <= b(p,ij)
                            
                            % Coefficient for h_{ik}
                            tb_row_index = [tb_row_index,countb];
                            tb_col_index = [tb_col_index,ik];
                            
                            coeff_hik = (1/2)*U(k,ij);
                            
                            tb_mat_val = [tb_col_index,kj];
                            
                            % Coefficient for h_{kj}
                            tb_row_index = [tb_row_index,countb];
                            tb_col_index = [tb_col_index,kj];
                            
                            coeff_hkj = (1/2)*((-1)^(p-1))*V(k,ij);
                            
                            tb_mat_val = [tb_mat_val,coeff_hkj];
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % Inequality: -t_p{ij} <= -a(p,ij)
                            
                            % Coefficient for h_{ik}
                            ta_row_index = [ta_row_index,counta];
                            ta_col_index = [ta_col_index,ik];
                            
                            coeff_hik = -(1/2)*U(k,ij);
                            
                            ta_mat_val = [ta_col_index,kj];
                            
                            % Coefficient for h_{kj}
                            ta_row_index = [ta_row_index,counta];
                            ta_col_index = [ta_col_index,kj];
                            
                            coeff_hkj = -(1/2)*((-1)^(p-1))*V(k,ij);
                            
                            ta_mat_val = [ta_mat_val,coeff_hkj];
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                        end % end of if/elseif statement, if(ik == kj)
                        
                    else
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Inequality: t_p{ij} <= b(p,ij)
                        
                        % For h_{kj}
                        kj = (k-1)*m + j; % k=1:n
                        
                        % k > m, only generate h_{kj}
                        tb_row_index = [tb_row_index,counta];
                        tb_col_index = [tb_col_index,kj];
                        
                        coeff_hkj = (1/2)*((-1)^(p-1))*V(k,ij);
                        
                        tb_mat_val = [tb_mat_val,coeff_hkj];
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Inequality: -t_p{ij} <= -a(p,ij)
                        
                        % For h_{kj}
                        kj = (k-1)*m + j; % k=1:n
                        
                        % k > m, only generate h_{kj}
                        ta_row_index = [ta_row_index,counta];
                        ta_col_index = [ta_col_index,kj];
                        
                        coeff_hkj = -(1/2)*((-1)^(p-1))*V(k,ij);
                        
                        ta_mat_val = [ta_mat_val,coeff_hkj];
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                    end
                    
                    if(k <= n && max(m,n) == m)
                        
                        % For h_{ik}
                        ik = (i-1)*m + k; % k=1:m
                        
                        % For h_{kj}
                        kj = (k-1)*m + j; % k=1:n
                        
                        if(ik == kj)
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % Inequality: t_p{ij} <= b(p,ij)
                            
                            % Coefficients for h_{ik} and h_{kj} combined
                            tb_row_index = [tb_row_index,countb];
                            tb_col_index = [tb_col_index,kj]; % or ik
                            
                            coeff_val = (1/2)*(U(k,ij) + ((-1)^(p-1))*V(k,ij))
                            
                            tb_mat_val = [tb_mat_val,coeff_val];
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % Inequality: -t_p{ij} <= -a(p,ij)
                            
                            % Coefficients for h_{ik} and h_{kj} combined
                            ta_row_index = [ta_row_index,counta];
                            ta_col_index = [ta_col_index,kj]; % or ik
                            
                            coeff_val = -(1/2)*(U(k,ij) + ((-1)^(p-1))*V(k,ij))
                            
                            ta_mat_val = [ta_mat_val,coeff_val];
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        elseif(ik ~= kj)
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % Inequality: t_p{ij} <= b(p,ij)
                            
                            % Coefficient for h_{ik}
                            tb_row_index = [tb_row_index,countb];
                            tb_col_index = [tb_col_index,ik];
                            
                            coeff_hik = (1/2)*U(k,ij);
                            
                            tb_mat_val = [tb_col_index,kj];
                            
                            % Coefficient for h_{kj}
                            tb_row_index = [tb_row_index,countb];
                            tb_col_index = [tb_col_index,kj];
                            
                            coeff_hkj = (1/2)*((-1)^(p-1))*V(k,ij);
                            
                            tb_mat_val = [tb_mat_val,coeff_hkj];
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                            % Inequality: -t_p{ij} <= -a(p,ij)
                            
                            % Coefficient for h_{ik}
                            ta_row_index = [ta_row_index,counta];
                            ta_col_index = [ta_col_index,ik];
                            
                            coeff_hik = (1/2)*U(k,ij);
                            
                            ta_mat_val = [ta_col_index,kj];
                            
                            % Coefficient for h_{kj}
                            ta_row_index = [ta_row_index,counta];
                            ta_col_index = [ta_col_index,kj];
                            
                            coeff_hkj = (1/2)*((-1)^(p-1))*V(k,ij);
                            
                            ta_mat_val = [ta_mat_val,coeff_hkj];
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                        end % end of if/elseif statement, if(ik == kj)
                        
                    else
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Inequality: t_p{ij} <= b(p,ij)
                        
                        % For h_{ik}
                        ik = (i-1)*m + k; % k=1:m
                        
                        % k > n, only generate h_{ik}
                        tb_row_index = [tb_row_index,countb];
                        tb_col_index = [tb_col_index,ik];
                        
                        coeff_hik = (1/2)*U(k,ij);
                        
                        tb_mat_val = [tb_mat_val,coeff_hik];
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Inequality: -t_p{ij} <= -a(p,ij)
                        
                        % For h_{ik}
                        ik = (i-1)*m + k; % k=1:m
                        
                        % k > n, only generate h_{ik}
                        ta_row_index = [ta_row_index,counta];
                        ta_col_index = [ta_col_index,ik];
                        
                        coeff_hik = -(1/2)*U(k,ij);
                        
                        ta_mat_val = [ta_mat_val,coeff_hik];
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                    end
                    
                end % end of for loop, for k=1:max(m,n)
                
            end % end of for loop, for j=1:m
        end % end of for loop, for i=1:n
        
    end % end of for loop, for p=1:2
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NOTE(about model.quadcon.[] in GUROBI):  x'(Qc)x + q'x <= b
    %
    %   Qc must be square, sparse; row/col. counts = # of cols. of A (the
    %   linear constraint matrix) : stored in model.quadcon.Qc
    %
    %   q vector defines linear terms of constraint; specify for each
    %   col. of A (dense vector) : stored in model.quadcon.q
    %
    %   b defines right-hand side : stored in model.quadcon.rhs
    
    % Quadratic Lifting Inequalities: h_{ij}, w_p{ij}, K_{ij}{kl} variables
    
    countql = 0;
    
    for p=1:2
        
        for i=1:n
            for j=1:m
                
                countql = countql + 1;
                
                ij = (i-1)*m + j; % Mapping index for H(i,j) - var. location
                
                pij = (mod(p,2) + 1)*n*m + ij; % Given p, index for w_(mod(p,2) + 1){ij} (w.r.t. t_p{ij})
                
                % p = 1 -> t_1{ij}, w_((p mod 2) + 1){ij} - w_2{ij}, K_{ij}{kl} positive
                % p = 2 -> t_2{ij}m w((p mod 2) + 1){ij} - w_1{ij}, K_{ij}{kl} negative
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Quadratic Term: h_{ik}h_{il}
                for k=1:m
                    for l=1:m
                        
                        ik = (i-1)*m + k; % index for h_{ik}
                        il = (i-1)*m + l; % Index for h_{il}
                        
                        if(ik == il)
                            
                            Qc_row_index = [Qc_row_index,ik]; % or il
                            Qc_col_index = [Qc_col_index,il]; % or ik
                            
                            coeff_hik_hil = (1/4)*U(k,ij)*U(l,ij);
                            
                            Qc_mat_val = [Qc_mat_val,coeff_hik_hil];
                            
                        elseif(ik ~= il)
                            
                            % Symmetric, split across h_{ik}h_{il} and h_{il}h_{ik}
                            
                            % h_{ik}h_{il} portion
                            Qc_row_index = [Qc_row_index,ik];
                            Qc_col_index = [Qc_col_index,il];
                            
                            coeff_hik_hil = (1/2)*((1/4)*U(k,ij)*U(l,ij));
                            
                            Qc_mat_val = [Qc_mat_val,coeff_hik_hil];
                            
                            % h_{il}h_{ik} portion
                            Qc_row_index = [Qc_row_index,il];
                            Qc_col_index = [Qc_col_index,ik];
                            
                            coeff_hil_hik = (1/2)*((1/4)*U(l,ij)*U(k,ij));
                            
                            Qc_mat_val = [Qc_mat_val,coeff_hil_hik];
                            
                        end % end of if/elseif statement, if(ik == il)
                        
                    end
                end
                
                % Quadratic Term: h_{ik}h_{lj}
                for k=1:m
                    for l=1:n
                        
                        ik = (i-1)*m + k; % index for h_{ik}
                        lj = (l-1)*m + j; % Index for h_{lj}
                        
                        if(ik == lj)
                            
                            Qc_row_index = [Qc_row_index,ik]; % or lj
                            Qc_col_index = [Qc_col_index,lj]; % or ik
                            
                            coeff_hik_hlj = (1/4)*((-1)^(p-1))*U(k,ij)*V(l,ij);
                            
                            Qc_mat_val = [Qc_mat_val,coeff_hik_hlj];
                            
                        elseif(ik ~= lj)
                            
                            % Symmetric, split across h_{ik}h_{il} and h_{il}h_{ik}
                            
                            % h_{ik}h_{il} portion
                            Qc_row_index = [Qc_row_index,ik];
                            Qc_col_index = [Qc_col_index,lj];
                            
                            coeff_hik_hlj = (1/2)*((-1)^(p-1))*((1/4)*U(k,ij)*V(l,ij));
                            
                            Qc_mat_val = [Qc_mat_val,coeff_hik_hlj];
                            
                            % h_{il}h_{ik} portion
                            Qc_row_index = [Qc_row_index,lj];
                            Qc_col_index = [Qc_col_index,ik];
                            
                            coeff_hlj_hik = (1/2)*((-1)^(p-1))*((1/4)*V(l,ij)*U(k,ij));
                            
                            Qc_mat_val = [Qc_mat_val,coeff_hlj_hik];
                            
                        end % end of if/elseif statement, if(ik == il)
                        
                    end
                end
                
                % Quadratic Term: h_{kj}h_{il}
                for k=1:n
                    for l=1:m
                        
                        kj = (k-1)*m + j; % index for h_{kj}
                        il = (i-1)*m + l; % Index for h_{il}
                        
                        if(kj == il)
                            
                            Qc_row_index = [Qc_row_index,kj]; % or il
                            Qc_col_index = [Qc_col_index,il]; % or kj
                            
                            coeff_hkj_hil = (1/4)*((-1)^(p-1))*V(k,ij)*U(l,ij);
                            
                            Qc_mat_val = [Qc_mat_val,coeff_hkj_hil];
                            
                        elseif(kj ~= il)
                            
                            % Symmetric, split across h_{ik}h_{il} and h_{il}h_{ik}
                            
                            % h_{ik}h_{il} portion
                            Qc_row_index = [Qc_row_index,kj];
                            Qc_col_index = [Qc_col_index,il];
                            
                            coeff_hkj_hil = (1/2)*((-1)^(p-1))*((1/4)*V(k,ij)*U(l,ij));
                            
                            Qc_mat_val = [Qc_mat_val,coeff_hkj_hil];
                            
                            % h_{il}h_{ik} portion
                            Qc_row_index = [Qc_row_index,il];
                            Qc_col_index = [Qc_col_index,kj];
                            
                            coeff_hil_hkj = (1/2)*((-1)^(p-1))*((1/4)*U(l,ij)*V(k,ij));
                            
                            Qc_mat_val = [Qc_mat_val,coeff_hil_hkj];
                            
                        end % end of if/elseif statement, if(ik == il)
                        
                    end
                end
                
                % Quadratic Term: h_{kj}h_{lj}
                for k=1:n
                    for l=1:n
                        
                        kj = (k-1)*m + j; % index for h_{kj}
                        lj = (l-1)*m + j; % Index for h_{lj}
                        
                        if(kj == lj)
                            
                            Qc_row_index = [Qc_row_index,kj]; % or lj
                            Qc_col_index = [Qc_col_index,lj]; % or kj
                            
                            coeff_hkj_hlj = (1/4)*V(k,ij)*V(l,ij);
                            
                            Qc_mat_val = [Qc_mat_val,coeff_hkj_hlj];
                            
                        elseif(kj ~= lj)
                            
                            % Symmetric, split across h_{ik}h_{il} and h_{il}h_{ik}
                            
                            % h_{ik}h_{il} portion
                            Qc_row_index = [Qc_row_index,kj];
                            Qc_col_index = [Qc_col_index,lj];
                            
                            coeff_hkj_hlj = (1/2)*((1/4)*V(k,ij)*V(l,ij));
                            
                            Qc_mat_val = [Qc_mat_val,coeff_hkj_hlj];
                            
                            % h_{il}h_{ik} portion
                            Qc_row_index = [Qc_row_index,lj];
                            Qc_col_index = [Qc_col_index,kj];
                            
                            coeff_hlj_hkj = (1/2)*((1/4)*V(l,ij)*V(k,ij));
                            
                            Qc_mat_val = [Qc_mat_val,coeff_hlj_hkj];
                            
                        end % end of if/elseif statement, if(ik == il)
                        
                    end
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Linear Term: K_{ij}{kl}
                
                % p = 1 -> linked to t_1{ij}, negative coefficient
                % p = 2 -> linked to t_2{ij}, positive coefficient
                
                % UPDATE(6/6/17): Incorporate Z (nm) into col. indices
                
                for k=1:m
                    for l=1:n
                        q_row_index = [q_row_index,countql];
                        
                        kl = (l-1)*m + k;
                        
                        if(ij <= kl)
                            
                            sum_t = sum(ones(1,1:ij-1));
                            
                            % Accounts for all Z(i,j), w_p{ij}, shift accordingly
                            q_col_index = [q_col_index,n*m + 2*n*m + (ij)*n*m + kl - sum_t];
                            
                        elseif(ij > kl)
                            
                            sum_t = sum(ones(1,1:kl-1));
                            
                            % Accounts for all Z(i,j), w_p{ij} [UPDATE: 6/6/17]
                            q_col_index = [q_col_index,n*m + 2*n*m + (kl)*n*m + ij - sum_t];
                            
                        end % end of if/elseif statement
                        
                        coeff_Kijkl = ((-1)^(p))*U(k,ij)*V(l,ij);
                        
                        q_mat_val = [q_mat_val,coeff_Kijkl];
                        
                    end
                end
                
                % Linear Term: w_p{ij}
                
                % p = 1 -> linked to t_1{ij}, so w_2{ij}
                % p = 2 -> linked to t_2{ij}, so w_1{ij}
                
                % NOTE: pij accounts for this relationship
                
                q_row_index = [q_row_index,countql];
                
                % Account for Z(i,j), shift accordingly [UPDATE: 6/6/17]
                q_col_index = [q_col_index,n*m + pij];
                q_mat_val = [q_mat_val,1];
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            end % end of for loop, for j=1:m
        end % end of for loop, for i=1:n
        
    end % end of for loop, for p=1:2
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Construct sparse representations
    
    %PT_sparse = sparse(PT_row_index,PT_col_index,PT_mat_val,(n*m)^2,2*n*m + nchoosek(n*m,2));
    
    % Secant constraints:
    
    
    % UPDATE(6/6/17): Incorporate Z (nm) variables (order: H,Z,w,K)
    secant_sparse = sparse(secant_row_index,secant_col_index,secant_mat_val,count_sec,n*m + n*m + 2*n*m + n*m + nchoosek(n*m,2));
    
    % Bounds on t_p{ij}:
    
    
    % UPDATE(6/6/17): Incorporate Z (nm) variables (order: H,Z,w,K)
    ta_sparse = sparse(ta_row_index,ta_col_index,ta_mat_val,counta,n*m + n*m + 2*n*m + n*m + nchoosek(n*m,2));
    tb_sparse = sparse(tb_row_index,tb_col_index,tb_mat_val,countb,n*m + n*m + 2*n*m + n*m + nchoosek(n*m,2));
    
    % Convex quadratics:
    
    % UPDATE(6/6/17): Incorporate Z (nm) variables (order: H,Z,w,K)
    Qc_sparse = sparse(Qc_row_index,Qc_col_index,Qc_mat_val,3*n*m + nchoosek(n*m,2),n*m + n*m + 2*n*m + n*m + nchoosek(n*m,2));
    
    % Linear Terms:
    
    % UPDATE(6/6/17): Incorporate Z (nm) variables (order: H,Z,w,K)
    q_sparse = sparse(q_row_index,q_col_index,q_mat_val,countql,n*m + n*m + 2*n*m + n*m + nchoosek(n*m,2));
    
else
    
    secant_sparse = [];
    
    ta_sparse = [];
    tb_sparse = [];
    
    Qc_sparse = [];
    q_sparse = [];
    
end % end of if statement, if(strcmp(set_prop(7),''))

% Box constraints on H: L(i,j) <= H(i,j) <= M(i,j)

if(strcmp(set_prop{8},'included'))
    
    for i=1:n
        for j=1:m
            
            ij = (i-1)*m + j; % Index for matrix entries (H,Mu,Ld)
            
            % Generate box constraints on H (using L/M)
            
            % Upper box constraint: H(i,j) <= Mu(i,j)
            %   - Constraint Type: rgl, '<'
            %   - rhs: Mu(i,j)
            
            Mu_row_index = [Mu_row_index,ij];
            Mu_col_index = [Mu_col_index,ij];
            Mu_mat_val = [Mu_mat_val,1];
            
            
            % Lower box constraint: -H(i,j) <= -Ld(i,j)
            %   - Constraint Type: rgl, '<'
            %   - rhs: -Ld(i,j)
            
            Ld_row_index = [Ld_row_index,ij];
            Ld_col_index = [Ld_col_index,ij];
            Ld_mat_val = [Ld_mat_val,-1];
            
            
        end % end of for loop, for j=1:m
    end % end of for loop, for i=1:n
    
    % Check for P1-tightening, Non-sym P2, McCormick Ineq, or Quadratic Ineq
    if(strcmp(set_prop{4},'included') || strcmp(set_prop{5},'included') || strcmp(set_prop{6},'included') || strcmp(set_prop{7},'included'))
        
        % Variable space: H,Z,w,K
        
        Mu_sparse = sparse(Mu_row_index,Mu_col_index,Mu_mat_val,n*m,(n*m + n*m + 2*n*m) + n*m + nchoosek(n*m,2));
        Ld_sparse = sparse(Ld_row_index,Ld_col_index,Ld_mat_val,n*m,(n*m + n*m + 2*n*m) + n*m + nchoosek(n*m,2));
        
    else
        
        % Variable space: H,Z
        
        Mu_sparse = sparse(Mu_row_index,Mu_col_index,Mu_mat_val,n*m,n*m + n*m);
        Ld_sparse = sparse(Ld_row_index,Ld_col_index,Ld_mat_val,n*m,n*m + n*m);
        
    end
    
else
    
    Mu_sparse = [];
    Ld_sparse = [];
    
end % end of if statement, if(strcmp(set_prop{8},''))


% Sparse representation - accumulate across all possible constraints
%keyboard
M_sparse = [one_NL_sparse; one_NU_sparse; P1_sparse; P3_sparse; P4_sparse; PT_sparse; NS_sparse; mccl_sparse; secant_sparse; tb_sparse; ta_sparse; Ld_sparse; Mu_sparse; P_sym_sparse;P2L_sparse];

% NOTE: Return Qc_sparse and q_sparse

end % end of function, constr_matrix_gen