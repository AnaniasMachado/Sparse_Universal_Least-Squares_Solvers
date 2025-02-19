clear;
  
rowop=[50 80 100] ;
colop=[50 80 100] ;
rankop=[0.10 0.50];
densop=[0.25 0.5 1.0];
numrep=5;

numrowcol=length(rowop);
numrank=length(rankop);
numden=length(densop);

M_aux = [];

for irowcol=1:numrowcol
    for irank = 1:numrank
        for iden = 1:numden
            for irep = 1:numrep   
                
                m=rowop(irowcol);
                n=colop(irowcol);
                r=rankop(irank)*n;
                dens=densop(iden);
                name_M = ['A' num2str(m) '_' num2str(n) '_' num2str(r) '_' num2str(dens) '_' num2str(irep) '.mat']
                load(name_M);
                
                A = full(A);
                
                name_H1 = ['H_P1_' num2str(m) '_' num2str(n) '_' num2str(r) '_' num2str(dens) '_' num2str(irep) '.mat'];
                name_H13 = ['H_P13_' num2str(m) '_' num2str(n) '_' num2str(r) '_' num2str(dens) '_' num2str(irep) '.mat'];
                name_H123 = ['H_P123_' num2str(m) '_' num2str(n) '_' num2str(r) '_' num2str(dens) '_' num2str(irep) '.mat'];
                name_H1sym = ['H_P1sym_' num2str(m) '_' num2str(n) '_' num2str(r) '_' num2str(dens) '_' num2str(irep) '.mat'];
                
                [norma_gurobi1,time1,H1] = gurobi_MainFunctionP1(A,m,n);
                [norma_gurobi13,time13,H13] = gurobi_MainFunctionP1P3(A,m,n);
                [norma_gurobi123,time123,H123] = gurobi_MainFunctionP1P2P3(A,m,n);
                [norma_gurobi1sym,time1sym,H1sym] = gurobi_MainFunctionP1sym(A,m,n);
                
                save(name_H1,'norma_gurobi1','time1','H1');
                save(name_H13,'norma_gurobi13','time13','H13');
                save(name_H123,'norma_gurobi123','time123','H123');
                save(name_H1sym,'norma_gurobi1sym','time1sym','H1sym');
                
                M_aux = [M_aux;[norma_gurobi1,norma_gurobi13,norma_gurobi123,norma_gurobi1sym,time1,time13,time123,time1sym]];
            end
        end
    end
end