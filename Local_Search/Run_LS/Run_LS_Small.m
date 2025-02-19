clear;
   
%Define the group of intances:

%small

rowop=[50 80 100] ;
colop=[50 80 100] ;
rankop=[0.10 0.50];
densop=[0.25 0.5 1.0];
numrep=5;


% %medium
% foldernamedata='medium/';
% if greedy == 1
%     foldernameresults='C:\Users\gabponte\Desktop\Testes Montreal\medium\Results_medium\Greedy\';
% else
%     foldernameresults='C:\Users\gabponte\Desktop\Testes Montreal\medium\Results_medium\Greedy_Light\';
% end
% rowop=[1000 2000] ;
% colop=[1000 2000] ;
% rankop=[0.05 0.10];
% densop=[0.25 0.5 1.0];
% numrep=30;

% %large
% foldernamedata='large/';
% if greedy == 1
%     foldernameresults='large/Results_large/Greedy/';
% else
%     foldernameresults='large/Results_lerge/HQR/';
% end
% rowop=[5000 10000] ;
% colop=[1000 1000] ;
% rankop=[0.05 0.10];
% densop=[1.0];
% numrep=3;

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
%                 name_M = ['A' num2str(m) '_' num2str(n) '_' num2str(r) '_' num2str(dens) '_' num2str(irep) '.mat']
%                 load(name_M)
                
                name_Hinit = [fileHinit num2str(m) '_' num2str(n) '_' num2str(r) '_' num2str(dens) '_' num2str(irep) '.mat'];
                load(name_Hinit);
                var = [var;[factor_mc,factor_mr]];
                %keyboard

%                 A=full(A);
%                 normH_Init=sum(sum(abs(inv(A(R,C)))));
%                 time_Init=time;
% 
%                 [normLSFI_Det,timeLSFI_Det,swaps_LSFI_Det,R_LSFI_Det,C_LSFI_Det] = LSFI_Det (A,r,m,n,R,C);
%                 name_H_LSFI_Det = ['H_LSFI_Det_' num2str(m) '_' num2str(n) '_' num2str(r) '_' num2str(dens) '_' num2str(irep) '.mat'];
%                 save([foldernameresults name_H_LSFI_Det],'R_LSFI_Det', 'C_LSFI_Det', 'swaps_LSFI_Det', 'timeLSFI_Det', 'normLSFI_Det');
% 
%                 [normLSLAFI_Det,timeLSLAFI_Det,swaps_LSLAFI_Det,R_LSLAFI_Det,C_LSLAFI_Det] = LSLAFI_Det (A,r,m,n,R,C);
%                 name_H_LSLAFI_Det = ['H_LSLAFI_Det_' num2str(m) '_' num2str(n) '_' num2str(r) '_' num2str(dens) '_' num2str(irep) '.mat'];
%                 save([foldernameresults name_H_LSLAFI_Det],'R_LSLAFI_Det', 'C_LSLAFI_Det', 'swaps_LSLAFI_Det', 'timeLSLAFI_Det', 'normLSLAFI_Det');
% 
%                 [normLSBI_Det,timeLSBI_Det,swaps_LSBI_Det,R_LSBI_Det,C_LSBI_Det] = LSBI_Det (A,r,m,n,R,C);
%                 name_H_LSBI_Det = ['H_LSBI_Det_' num2str(m) '_' num2str(n) '_' num2str(r) '_' num2str(dens) '_' num2str(irep) '.mat'];
%                 save([foldernameresults name_H_LSBI_Det],'R_LSBI_Det', 'C_LSBI_Det', 'swaps_LSBI_Det', 'timeLSBI_Det', 'normLSBI_Det');
% 
%                 [normLSFI_Det_FI_Norm,timeLSFI_Det_FI_Norm,swaps_LSFI_Det_FI_Norm,R_LSFI_Det_FI_Norm,C_LSFI_Det_FI_Norm] = LSFI_Norm(A,r,m,n,R_LSFI_Det,C_LSFI_Det);
%                 name_H_LSFI_Det_FI_Norm = ['H_LSFI_Det_FI_Norm_' num2str(m) '_' num2str(n) '_' num2str(r) '_' num2str(dens) '_' num2str(irep) '.mat'];
%                 save([foldernameresults name_H_LSFI_Det_FI_Norm],'R_LSFI_Det_FI_Norm', 'C_LSFI_Det_FI_Norm', 'swaps_LSFI_Det_FI_Norm', 'timeLSFI_Det_FI_Norm', 'normLSFI_Det_FI_Norm');
% 
%                 [normLSLAFI_Det_FI_Norm,timeLSLAFI_Det_FI_Norm,swaps_LSLAFI_Det_FI_Norm,R_LSLAFI_Det_FI_Norm,C_LSLAFI_Det_FI_Norm] = LSFI_Norm(A,r,m,n,R_LSLAFI_Det,C_LSLAFI_Det);
%                 name_H_LSLAFI_Det_FI_Norm = ['H_LSLAFI_Det_FI_Norm_' num2str(m) '_' num2str(n) '_' num2str(r) '_' num2str(dens) '_' num2str(irep) '.mat'];
%                 save([foldernameresults name_H_LSLAFI_Det_FI_Norm],'R_LSLAFI_Det_FI_Norm', 'C_LSLAFI_Det_FI_Norm', 'swaps_LSLAFI_Det_FI_Norm', 'timeLSLAFI_Det_FI_Norm', 'normLSLAFI_Det_FI_Norm');
% 
%                 [normLSBI_Det_FI_Norm,timeLSBI_Det_FI_Norm,swaps_LSBI_Det_FI_Norm,R_LSBI_Det_FI_Norm,C_LSBI_Det_FI_Norm] = LSFI_Norm(A,r,m,n,R_LSBI_Det,C_LSBI_Det);
%                 name_H_LSBI_Det_FI_Norm = ['H_LSBI_Det_FI_Norm_' num2str(m) '_' num2str(n) '_' num2str(r) '_' num2str(dens) '_' num2str(irep) '.mat'];
%                 save([foldernameresults name_H_LSBI_Det_FI_Norm],'R_LSBI_Det_FI_Norm', 'C_LSBI_Det_FI_Norm', 'swaps_LSBI_Det_FI_Norm', 'timeLSBI_Det_FI_Norm', 'normLSBI_Det_FI_Norm');
% 
%                 if rowtable == 1
%                     var=[m,n,r,dens,irep,normH_Init,normLSFI_Det,normLSLAFI_Det,normLSBI_Det,normLSFI_Det_FI_Norm,normLSLAFI_Det_FI_Norm,normLSBI_Det_FI_Norm,time_Init,timeLSFI_Det,timeLSLAFI_Det,timeLSBI_Det,timeLSFI_Det_FI_Norm,timeLSLAFI_Det_FI_Norm,timeLSBI_Det_FI_Norm,swaps_LSFI_Det,swaps_LSLAFI_Det,swaps_LSBI_Det,swaps_LSFI_Det_FI_Norm, swaps_LSLAFI_Det_FI_Norm,swaps_LSBI_Det_FI_Norm];
%                 else
%                     var=[var;m,n,r,dens,irep,normH_Init,normLSFI_Det,normLSLAFI_Det,normLSBI_Det,normLSFI_Det_FI_Norm,normLSLAFI_Det_FI_Norm,normLSBI_Det_FI_Norm,time_Init,timeLSFI_Det,timeLSLAFI_Det,timeLSBI_Det,timeLSFI_Det_FI_Norm,timeLSLAFI_Det_FI_Norm,timeLSBI_Det_FI_Norm,swaps_LSFI_Det,swaps_LSLAFI_Det,swaps_LSBI_Det,swaps_LSFI_Det_FI_Norm, swaps_LSLAFI_Det_FI_Norm,swaps_LSBI_Det_FI_Norm];
%                 end
%                 rowtable=rowtable+1;
            end
        end
    end
end

folder = 'Results';
if ~exist(folder, 'dir')
    mkdir(folder);
end
fullFileName = fullfile(folder, baseFileName);
text={'m','n','r','dens','inst','|H_Init|','|FI_Det|','|LAFI_Det|','|BI_Det|','|FI_Det_FI_Norm|','|LAFI_Det_FI_Norm|','|BI_Det_FI_Norm|','t_Init','t_FI_Det','t_LAFI_Det','t_LSBI_Det','t_FI_Det_FI_Norm','t_LAFI_Det_FI_Norm','t_BI_Det_FI_Norm','swaps_FI_Det','swaps_LAFI_Det','swaps_BI_Det','swaps_FI_Det_FI_Norm', 'swaps_LAFI_Det_FI_Norm','swaps_BI_Det_FI_Norm'};
xlswrite(fullFileName, text,1,'A1');
xlswrite(fullFileName ,var,1,'A2');


    