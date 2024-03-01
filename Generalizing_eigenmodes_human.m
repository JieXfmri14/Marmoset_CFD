function [recon_ratio,recon_FC_r]= Generalizing_eigenmodes_human(mypath,n_null)
% The CC matrix from the homologous brain regions was used to reconstruct the human Bold-fMRI

%% load homology CC matrix and HCP Bold-fMRI
% homology CC matrix
W=load(strcat(mypath,'homology_marmoset_human/homology_matrix_11ROI.mat'));      
W=W.homology_matrix; 
n_ROI = size(W,2);  

%% HCP Bold-fMRI, left brain
X_RS=load(strcat(mypath,'homology_marmoset_human/HCP_BOLD_demo.mat'));
X_RS = X_RS.MMPsignal;
n_sub=size(X_RS,3);

% Extract homology signals
MMP = xlsread(strcat(mypath,'homology_marmoset_human/Human_MMP_homology_label_11ROI.xlsx'),'B2:C181');
raw_ROI_label = MMP(:,1);
MMP_homology = MMP(:,2);    
order = get_homology_MMP(raw_ROI_label,MMP_homology);
X_RS_homology = X_RS(order,:,:);
% dmean mean centering
mean_data = mean(X_RS_homology,2);
zX_RS = X_RS_homology - mean_data;       
RX_RS = X_RS_homology; 

% non-homology singals
nonhomology = raw_ROI_label;
nonhomology(order) = [];
X_RS_nonhomology = X_RS(nonhomology,:,:);
% dmean mean centering
mean_data = mean(X_RS_nonhomology,2);
zX_RS_nonhomology = X_RS_nonhomology - mean_data;     
RX_RS_nonhomology = X_RS_nonhomology;                           

%% Laplacian Decomposition
L = Computer_laplacian_matrix(W);
[U_old,Lambda] = eig(L);   
[Lambda, Ind]=sort(diag(Lambda));     
U=U_old(:,Ind);

% =========================================================================
%     Calculate reconstruction accuracy using 1 to num_modes eigenmodes    
% =========================================================================
time_point = size(zX_RS,2);
n_modes = n_ROI;
X_all = zeros(n_modes,n_ROI,time_point,n_sub);
N_all = zeros(n_modes,n_ROI,n_sub);
%% reconstruction brain activity(BOLD-fMRI)
for mode= 1:n_modes
    for s=1:n_sub
        % Calculate reconstruction beta coefficients
        X_hat(:,:,s)=U'*zX_RS(:,:,s); 
        % reconstructed BOLD-fMRI
        M=zeros(size(U));
        M(:,1:mode)=U(:,1:mode);   
        X_all(mode,:,:,s)=M*X_hat(:,:,s);  
        
       %% norms of reconstructed BOLD-fMRI
        for r=1:n_ROI
            temp = X_all(mode,r,:,s);
            recon_signal_per = reshape(temp,time_point,1);
            N_all(mode,r,s) = norm(recon_signal_per);   
            
            %total energy
            signal(r,s)  = norm(zX_RS(r,:,s));  
        end
    end
end
recon_signal = mean(N_all,3);                      
recon_energy = mean(recon_signal,2);                       

acooss_sub = mean(signal,2);                      
real_energy = mean(acooss_sub);                        
recon_ratio = recon_energy/real_energy;

%% reconstruction FC
X_all = zeros(n_modes,n_ROI,time_point,n_sub);
for mode= 1:n_modes
    for s=1:n_sub
        X_hat(:,:,s)=U'*RX_RS(:,:,s); 
        M=zeros(size(U));
        M(:,1:mode)=U(:,1:mode);   
        X_all(mode,:,:,s)=M*X_hat(:,:,s);  
        
       %% FC of empirical signals
        FCpacereal(:,:,s)=corr(RX_RS(:,:,s)');
        clear temp
        temp(:,:,s) = X_all(mode,:,:,s);
        recon_FC(mode,:,:,s)=corr(temp(:,:,s)'); 
    end
    %Reconstruction FC
    group_FCpacereal = mean(FCpacereal,3); 
    group_recon_FCpacereal(mode,:,:) = mean(recon_FC(mode,:,:,:),4); 
    
    ind=find(triu(ones(n_ROI,n_ROI),+1)==1); 
    t1=squeeze(group_FCpacereal);
    t2=squeeze(group_recon_FCpacereal(mode,:,:));
    [recon_FC_r(mode,1),recon_FC_p(mode,1)]=corr(t1(ind),t2(ind));    
end
 
output  = [(1:n_modes)',recon_ratio,recon_FC_r];
xlswrite('Results/generalizing_CC_human/Human_homology_CC_recon_Accuracy.xlsx',output,strcat('A2:C',num2str(n_modes+1)))

%% null model
% =========================================================================
% (1)  Generate 1000 random digraphs (degree-preserving surrogate connectoms) further calculate reconstruction accuracy
% =========================================================================
[rewired_ratio, rewired_FC_r] = homology_CC_rewired_digraphs_Human(W,zX_RS,RX_RS,n_null)

% save data for R visual
for i = 1: n_modes
    data(1+n_null*(i-1):n_null*i,1) = i;   
    data(1+n_null*(i-1):n_null*i,2) = rewired_ratio(:,i);  
    data(1+n_null*(i-1):n_null*i,3) = rewired_FC_r(:,i);  
end
xlswrite('Results/generalizing_CC_human/Human_CC_recon_null_rewired.xlsx',data,strcat('A2:C',num2str(n_modes*n_null+1)))

% =========================================================================
% (2)  Generate 1000 non-homology regions' signals
% =========================================================================
[null_hom_ratio, null_hom_FC_r] = homology_non_signals_Human(U,zX_RS_nonhomology,RX_RS_nonhomology,n_null)

% save data for R visual
for i = 1: n_modes
    data(1+n_null*(i-1):n_null*i,1) = i;   
    data(1+n_null*(i-1):n_null*i,2) = null_hom_ratio(:,i);  
    data(1+n_null*(i-1):n_null*i,3) = null_hom_FC_r(:,i);  
end
xlswrite('Results/generalizing_CC_human/Human_CC_recon_null_nonhomology.xlsx',data,strcat('A2:C',num2str(n_modes*n_null+1)))

% =========================================================================
% (3)  Calculating p-value
% =========================================================================
for i =1:n_modes
    P_rewired_ratio(i) = sum(rewired_ratio(:,i)>recon_ratio(i))/n_null;
    P_rewired_FC(i) = sum(rewired_FC_r(:,i)> recon_FC_r(i))/n_null;
    
    P_nohomology_ratio(i) = sum(null_hom_ratio(:,i)>recon_ratio(i))/n_null;
    P_nohomology_FC(i) = sum(null_hom_FC_r(:,i)> recon_FC_r(i))/ n_null;
end
P_rewired_ratio(n_modes) = 0; P_rewired_FC(n_modes) = 0;
P_nohomology_ratio(n_modes) = 0; P_nohomology_FC(n_modes) = 0;

null_P = [P_rewired_ratio;P_rewired_FC;P_nohomology_ratio;P_nohomology_FC]';
xlswrite('Results/generalizing_CC_human/Human_CC_recon_null_Pvalue.xlsx',null_P)

