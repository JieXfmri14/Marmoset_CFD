function [N_low,N_high,mean_low,mean_high] = GSanalysis(zX_RS,Vhigh,Vlow,U)
% Graph signal analysis

data=zX_RS;             
n_ROI = size(data,1);
nsubjs_RS = size(data,3);

%% compute filter low- and high-frequency components
for s=1:nsubjs_RS 
    % Calculate reconstruction beta coefficients
    X_hat(:,:,s)=U'*data(:,:,s); 
    % reconstructed BOLD-fMRI
    X_low(:,:,s)=Vlow*X_hat(:,:,s);     
    X_high(:,:,s)=Vhigh*X_hat(:,:,s);       
    
    %% norms of the filtered signals across time points
    for r=1:n_ROI
        N_low(r,s)=norm(X_low(r,:,s));    
        N_high(r,s)=norm(X_high(r,:,s));
    end
end
%% mean across subjects
mean_low=mean(N_low,2); 
mean_high=mean(N_high,2); 

xlswrite('Results/filtered_low_high_signals/filtered_low_high_components.xlsx',[mean_low,mean_high],strcat('A2:B',num2str(n_ROI+1)))


end
