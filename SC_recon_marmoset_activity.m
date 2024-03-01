function [SC_recon_ratio, SC_recon_FC_r,...
    SC_rewired_ratio,SC_rewired_FC_r,...
    SC_moran_ratio,SC_moran_FC_r]= SC_recon_marmoset_activity(zX_RS,n_ROI,n_null)

SC = load('Data/marmoset_invivo_SC_Order.mat');
SC = SC.matrix_sc;
% only analysis 55 regions
SC = SC(1:n_ROI,1:n_ROI);

%% Threshold SC matrix to matched the density of CC 
% CC matrix: K = 1854; Kden = 0.6242
K = 1854;
Kden = 0.6242;
SC_com = reshape(SC,n_ROI*n_ROI,1);
SC_com = sort(SC_com,'descend');
thr = SC_com(K);
SC(SC<thr) = 0;
[kden,N,K] = density_dir(SC)

% =========================================================================
%  Calculate reconstruction accuracy using 1 to num_modes eigenmodes    
% =========================================================================
nsubjs_RS=size(zX_RS,3);
time_point = size(zX_RS,2);
n_modes = n_ROI;
X_all = zeros(n_modes,n_ROI,time_point,nsubjs_RS);
N_all = zeros(n_modes,n_ROI,nsubjs_RS);

%% Normalized graph Laplacian of SC 
L = Computer_laplacian_matrix(SC);

%% Laplacian Decomposition
[U_old,Lambda] = eig(L);                
[Lambda, Ind]=sort(diag(Lambda));     %lambda1<lambda2<lambda3
U=U_old(:,Ind);
%Verify orthogonality
orth_U = U*U';

%% reconstruction brain activity(BOLD-fMRI)
for mode= 1:n_modes
    for s=1:nsubjs_RS
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
            N_all(mode,r,s)=norm(recon_signal_per);   
            signal(r,s)  = norm(zX_RS(r,:,s));  
        end
    end
end
recon_signal = mean(N_all,3);                       
recon_signal_all = mean(recon_signal,2);            

acooss_sub = mean(signal,2);    
real_energy = mean(acooss_sub);  
SC_recon_ratio = recon_signal_all/real_energy;

%% reconstruction FC
for mode= 1:n_modes
    for s=1:nsubjs_RS       
       %% FC of empirical signals
        FCpacereal(:,:,s)=corr(zX_RS(:,:,s)');
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
    [SC_recon_FC_r(mode,1),SC_recon_FC_p(mode,1)]=corr(t1(ind),t2(ind));   
end

% % =========================================================================
% % Generate 1000 random digraphs (degree-preserving surrogate connectoms) further calculate reconstruction accuracy
% % =========================================================================
[SC_rewired_ratio, SC_rewired_FC_r]= SC_rewired_undigraphs(SC,zX_RS,n_null);

% % =========================================================================
% % Generate 1000 random Bold-signal further calculate reconstruction accuracy
% % =========================================================================
[SC_moran_ratio, SC_moran_FC_r]= CC_moran_signals(U,zX_RS,n_null);
