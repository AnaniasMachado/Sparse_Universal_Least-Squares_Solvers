clear;
  
rowop=[50 80 100] ;
colop=[50 80 100] ;
rankop=[0.10 0.50];
densop=[0.25 0.5 1.0];
numrep=5;

numrowcol=length(rowop);
numrank=length(rankop);
numden=length(densop);


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
                
                [R,C,time] = nsub(A,r);
                
                name_Hinit = ['H_nsub_' num2str(m) '_' num2str(n) '_' num2str(r) '_' num2str(dens) '_' num2str(irep) '.mat'];
                
                
                save(name_Hinit,'R','C','time');
            end
        end
    end
end