function [null_signal_energy, null_signal_recon_FC_r]= CC_permuted_signals(U,signals,null_graph)

% =========================================================================
%  Generate 1000 permuted functional activity by randomly permuting 
%  the raw rs-fMRI time series across brain regions further calculate reconstruction accuracy
% =========================================================================
n_ROI = size(U,1);
num_modes = n_ROI;
time_point = size(signals,2);
nsubjs_RS = size(signals,3);

X_all = zeros(num_modes,n_ROI,time_point,nsubjs_RS);
N_all = zeros(num_modes,n_ROI,nsubjs_RS);
null_signal_energy = zeros(null_graph,num_modes);
null_signal_recon_FC_r = zeros(null_graph,num_modes);
recon_FC_p = zeros(null_graph,num_modes);
null_ID = zeros(null_graph,n_ROI);

for null = 1:null_graph
    null
    permID = randperm(n_ROI,n_ROI);
    zX_RS = signals(permID',:,:);
    null_ID(null,:)  = permID';   

    for mode= 1:num_modes
        for s=1:nsubjs_RS        
            X_hat(:,:,s)=U'*zX_RS(:,:,s); 
           %% recon activity
            M=zeros(size(U));
            M(:,1:mode)=U(:,1:mode);   
            X_all(mode,:,:,s)=M*X_hat(:,:,s);  
            
           %% FC of empirical signals
            FCpacereal(:,:,s)=corr(zX_RS(:,:,s)');
            clear temp
            temp(:,:,s) = X_all(mode,:,:,s);
            recon_FC(mode,:,:,s)=corr(temp(:,:,s)'); 

           %% norms of reconstructed BOLD-fMRI
            for r=1:n_ROI
                temp = X_all(mode,r,:,s);
                recon_signal_per = reshape(temp,time_point,1);
                N_all(mode,r,s)=norm(recon_signal_per);   
                signal(r,s)  = norm(zX_RS(r,:,s));   
            end
        end
        
        %% Recon FC
        group_FCpacereal = mean(FCpacereal,3); 
        group_recon_FCpacereal(mode,:,:) = mean(recon_FC(mode,:,:,:),4); 

        ind=find(triu(ones(n_ROI,n_ROI),+1)==1); 
        t1=squeeze(group_FCpacereal);
        t2=squeeze(group_recon_FCpacereal(mode,:,:));
        [null_signal_recon_FC_r(null,mode),recon_FC_p(null,mode)]=corr(t1(ind),t2(ind));    
    end
    recon_signal_null = mean(N_all,3);                           
    recon_signal_all_null = mean(recon_signal_null,2);            
    
    acooss_sub = mean(signal,2);   
    real_energy = mean(acooss_sub);  
        
    null_signal_energy(null,:) = recon_signal_all_null/real_energy;  
end


