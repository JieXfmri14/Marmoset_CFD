function [W,U]= GSP_Laplacian(mypath)

%% load CC matrix
W=load(strcat(mypath,'Data/marmoset_CC_connectivity.mat')); 
W=W.matrix_fle; 

% only analysis 55 regions, source (55 regions) * target (55 regions)
n_ROI = size(W,2);  
W=W(1:n_ROI,:);
[kden,N,K] = density_dir(W);

% % =========================================================================
% %   Normalized directed graph Laplacian of the CC 
% % =========================================================================
L = Computer_laplacian_matrix(W);

% % =========================================================================
% %   Laplacian Decomposition
% % =========================================================================
[U_old,LambdaL] = eig(L);               
[LambdaL, IndL]=sort(diag(LambdaL));     %lambda1<lambda2<lambda3
U = U_old(:,IndL);
%Verify orthogonality
orth_U = U*U';

%% Frequency interpretation
%% Compute weighted zero crossings for Laplacian eigenvectors
for u=1:n_ROI              
    UU=U(:,u);             
    summ=0;
    for i=1:n_ROI-1          
        for j=i+1:n_ROI
            if (UU(i)*UU(j))<0                    
                summ=summ+(W(i,j)>0.01);           
            end
            wZC(u)=summ;         % Compute weighted zero of the Laplace eigenvectors
        end
    end
end

end
