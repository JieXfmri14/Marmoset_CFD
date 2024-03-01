function [recon_energy,recon_FC_r]= homology_CC_rewired_digraphs_Human(W,zX_RS,RX_RS,n_null)

% =========================================================================
%  Generate 1000 random digraphs (degree-preserving surrogate connectoms) further calculate reconstruction accuracy
% =========================================================================
n_ROI = size(W,1);
num_modes = n_ROI;
time_point = size(zX_RS,2);
nsubjs_RS = size(zX_RS,3);

X_all = zeros(num_modes,n_ROI,time_point,nsubjs_RS);
N_all = zeros(num_modes,n_ROI,nsubjs_RS);
recon_energy = zeros(n_null,num_modes);
recon_FC_r = zeros(n_null,num_modes);
recon_FC_p = zeros(n_null,num_modes);

for null = 1:n_null
    null
    [W_null, R] = null_model_dir_sign(W);   
    L = Computer_laplacian_matrix(W_null);
    %% Laplacian Decomposition
    [U_old,LambdaL] = eig(L);   
    [LambdaL, IndL]=sort(diag(LambdaL));     
    U=U_old(:,IndL);
    
    %% recon activity
    for mode= 1:num_modes
        for s=1:nsubjs_RS                  
            X_hat(:,:,s)=U'*zX_RS(:,:,s);      
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
    recon_signal_null = mean(N_all,3);                            
    recon_signal_all_null = mean(recon_signal_null,2);            
    
    acooss_sub = mean(signal,2);    
    real_energy = mean(acooss_sub);  
    recon_energy(null,:) = recon_signal_all_null/real_energy;
    
    %% recon FC
    for mode= 1:num_modes
        for s=1:nsubjs_RS                  
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
        % Reconstruction FC
        group_FCpacereal = mean(FCpacereal,3); 
        group_recon_FCpacereal(mode,:,:) = mean(recon_FC(mode,:,:,:),4); 
        
        ind=find(triu(ones(n_ROI,n_ROI),+1)==1); 
        t1=squeeze(group_FCpacereal);
        t2=squeeze(group_recon_FCpacereal(mode,:,:));
        [recon_FC_r(null,mode),recon_FC_p(null,mode)]=corr(t1(ind),t2(ind));    
    end 
end




