function [recon_energy, recon_FC_r]= CC_moran_signals(U,signals,n_null)

% =========================================================================
%  Generate 1000 permuted functional activity by Moran_randomization
%  the raw rs-fMRI time series across brain regions further calculate reconstruction accuracy
% =========================================================================
n_ROI = size(U,1);
n_modes = n_ROI;
time_point = size(signals,2);
nsubjs_RS = size(signals,3);

% X_all = zeros(num_modes,n_ROI,time_point,nsubjs_RS);
N_all = zeros(n_modes,n_ROI,nsubjs_RS);
recon_energy = zeros(n_null,n_modes);
recon_FC_r = zeros(n_null,n_modes);
recon_FC_p = zeros(n_null,n_modes);
null_ID = zeros(n_null,n_ROI);

%A matrix denoting distance between features
distance = load('Data/marmoset_distance_matrix.mat');
MEM = compute_mem(distance.matrix_distance55);

signal_MSR = zeros(n_ROI,time_point,n_null,nsubjs_RS);
for s=1:nsubjs_RS 
    signal_MSR(:,:,:,s) = moran_randomization(signals(:,:,s),MEM,n_null, ...
    'procedure','singleton','joint',true,'random_state',0);
end

for null = 1:n_null
    null
    
    zX_RS = signal_MSR(:,:,null,:);
    zX_RS = squeeze(zX_RS(:,:,1,:));
       
    for mode= 1:n_modes
        for s=1:nsubjs_RS  
            MSR_RS = zX_RS(:,:,s);
            MSR_RS(:, any(isnan(MSR_RS))) = [];
           
            X_hat=U'*MSR_RS; 
           %% recon activity
            M=zeros(size(U));
            M(:,1:mode)=U(:,1:mode);   
            X_all=M*X_hat;  
            
           %% FC of empirical signals
            FCpacereal(:,:,s)=corr(MSR_RS');
            clear temp
            recon_FC(mode,:,:,s)=corr(X_all'); 

           %% norms of reconstructed BOLD-fMRI
           for r=1:n_ROI
                N_all(mode,r,s)=norm(X_all(r,:));   
                signal(r,s)  = norm(MSR_RS(r,:));   
            end
        end
        
        %% Recon FC
        group_FCpacereal = mean(FCpacereal,3); 
        group_recon_FCpacereal(mode,:,:) = mean(recon_FC(mode,:,:,:),4); 

        ind=find(triu(ones(n_ROI,n_ROI),+1)==1); 
        t1=squeeze(group_FCpacereal);
        t2=squeeze(group_recon_FCpacereal(mode,:,:));
        [recon_FC_r(null,mode),recon_FC_p(null,mode)]=corr(t1(ind),t2(ind));    
    end
    recon_signal_null = mean(N_all,3);                           
    recon_signal_all_null = mean(recon_signal_null,2);            
    
    acooss_sub = mean(signal,2);   
    real_energy = mean(acooss_sub);  
        
    recon_energy(null,:) = recon_signal_all_null/real_energy;
end


