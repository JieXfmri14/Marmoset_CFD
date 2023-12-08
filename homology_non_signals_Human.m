function [null_hom_ratio, null_hom_recon_FC_r]= homology_non_signals_Human(U,demean_signals,RX_signals,null_graph)

% % =========================================================================
% %  Generate 1000 random non-homology regions' signal further calculate reconstruction accuracy
% % =========================================================================
n_ROI = size(U,1);
num_modes = n_ROI;
non_ROI = size(demean_signals,1);
time_point = size(demean_signals,2);
nsubjs_RS = size(demean_signals,3);

X_all = zeros(num_modes,n_ROI,time_point,nsubjs_RS);
N_all = zeros(num_modes,n_ROI,nsubjs_RS);
null_hom_recon_FC_r = zeros(null_graph,num_modes);
recon_FC_p = zeros(null_graph,num_modes);
null_ID = zeros(null_graph,n_ROI);

for null = 1:null_graph
    null
    % Generate 1000 random non-homology signal
    permID = randperm(non_ROI,n_ROI);
    null_ID(null,:)  = permID;    
    zX_RS = demean_signals(permID,:,:);     
    RX_RS = RX_signals(permID,:,:);      

    
    %% energy
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
        
    null_hom_ratio(null,:) = recon_signal_all_null/real_energy;
    
   %% recon FC
    for mode= 1:num_modes
        for s=1:nsubjs_RS                  
            X_hat(:,:,s)=U'*RX_RS(:,:,s);      
            M=zeros(size(U));
            M(:,1:mode)=U(:,1:mode);   
           %% reconstruct back full signal
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

        ind=find(triu(ones(n_ROI,n_ROI),+1)==1); % find indices of all values above the diagonal
        t1=squeeze(group_FCpacereal);
        t2=squeeze(group_recon_FCpacereal(mode,:,:));
        [null_hom_recon_FC_r(null,mode),recon_FC_p(null,mode)]=corr(t1(ind),t2(ind));    % calc r
    end
end
