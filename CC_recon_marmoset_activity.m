function [CC_recon_ratio,CC_recon_FC_r]= CC_recon_marmoset_activity(W,X_RS,U,n_null)
% (1) assess CC eigenmodes explain brain activity in marmosets
% (2) compare the accuracy of CC eigenmodes against SC eigenmodes

zX_RS = X_RS;

% % =========================================================================
% %  (1) Calculate reconstruction accuracy using 1 to num_modes eigenmodes    
% % =========================================================================
time_point = size(zX_RS,2);
zX_RS = zX_RS(:,1:time_point,:);
n_ROI = size(W,1);
n_modes = n_ROI;
nsubjs_RS = size(zX_RS,3);
X_all = zeros(n_modes,n_ROI,time_point,nsubjs_RS);
N_all = zeros(n_modes,n_ROI,nsubjs_RS);
%% reconstruction brain activity(BOLD-fMRI)
for mode= 1:n_modes
    for s=1:nsubjs_RS
        % Calculate reconstruction beta coefficients
        X_hat(:,:,s)=U'*zX_RS(:,:,s); 
        % reconstructed BOLD-fMRI
        M=zeros(size(U));
        M(:,1:mode)=U(:,1:mode);   
        X_all(mode,:,:,s)=M*X_hat(:,:,s);  
         
       %% norms (energy) of reconstructed BOLD-fMRI
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
CC_recon_ratio = recon_signal_all/real_energy;

%% reconstruction FC
for mode= 1:n_modes
    for s=1:nsubjs_RS       
       %% FC of empirical signals
        FCpacereal(:,:,s)=corr(zX_RS(:,:,s)');
        clear temp
        temp(:,:,s) = X_all(mode,:,:,s);
        recon_FC(mode,:,:,s)=corr(temp(:,:,s)'); 
    end
    % Reconstruction FC
    group_FCpacereal = mean(FCpacereal,3); 
    group_recon_FCpacereal(mode,:,:) = mean(recon_FC(mode,:,:,:),4); 
    
    % find indices of all values above the diagonal
    ind=find(triu(ones(n_ROI,n_ROI),+1)==1); 
    t1=squeeze(group_FCpacereal);
    t2=squeeze(group_recon_FCpacereal(mode,:,:));
    [CC_recon_FC_r(mode,1),CC_recon_FC_p(mode,1)]=corr(t1(ind),t2(ind)); 
end
CC_output  = [(1:n_modes)',CC_recon_ratio,CC_recon_FC_r];
xlswrite('Results/CC_recon_marmoset_activity/Marmoset_CC_recon_Accuracy.xlsx',CC_output,strcat('A2:C',num2str(n_modes+1)))

for mode= 1:n_modes
    recon_ratio_AUC(mode,1)=trapz(CC_recon_ratio(1:mode,1))/mode;
    recon_FC_r_AUC(mode,1)=trapz(CC_recon_FC_r(1:mode,1))/mode;
end
CC_outputAUC  = [(1:n_modes)',recon_ratio_AUC,recon_FC_r_AUC];
xlswrite('Results/CC_recon_marmoset_activity/Marmoset_CC_recon_Accuracy_AUC.xlsx',CC_outputAUC,strcat('A2:C',num2str(n_modes+1)))

%% null model
% % =========================================================================
% %  1) Generate 1000 random digraphs (degree-preserving surrogate connectomes) further calculate reconstruction accuracy
% % =========================================================================
[CC_rewired_ratio, CC_rewired_FC_r]= CC_rewired_digraphs(W,zX_RS,n_null);

% save data for R visual
for i = 1: n_modes
    data(1+n_null*(i-1):n_null*i,1) = i;   
    data(1+n_null*(i-1):n_null*i,2) = CC_rewired_ratio(:,i);  
    data(1+n_null*(i-1):n_null*i,3) = CC_rewired_FC_r(:,i);  
end
xlswrite('Results/CC_recon_marmoset_activity/Marmoset_CC_recon_null_rewired.xlsx',data,strcat('A2:C',num2str(n_modes*n_null+1)))

% % =========================================================================
% % 2) Generate 1000 Moran surrogate signals further calculate reconstruction accuracy
% % =========================================================================
[CC_moran_ratio, CC_moran_FC_r]= CC_moran_signals(U,zX_RS,n_null);
% save data for R visual
for i = 1: n_modes
    data(1+n_null*(i-1):n_null*i,1) = i;   
    data(1+n_null*(i-1):n_null*i,2) = CC_moran_ratio(:,i);  
    data(1+n_null*(i-1):n_null*i,3) = CC_moran_FC_r(:,i);  
end
xlswrite('Results/CC_recon_marmoset_activity/Marmoset_CC_recon_null_Moran.xlsx',data,strcat('A2:C',num2str(n_modes*n_null+1)))

% %=========================================================================================
% %  (2) Compare the performance of marmoset CC and invivo dMRI-based SC
% % ========================================================================================
[SC_recon_ratio, SC_recon_FC_r,...
    SC_rewired_ratio,SC_rewired_FC_r,...
    SC_moran_ratio,SC_moran_FC_r]= SC_recon_marmoset_activity(zX_RS,n_ROI,n_null)
output  = [(1:n_modes)',SC_recon_ratio,SC_recon_FC_r];
xlswrite('Results/CC_recon_marmoset_activity/Marmoset_SC_recon_Accuracy.xlsx',output,strcat('A2:C',num2str(n_modes+1)))

for mode= 1:n_modes
    SC_recon_ratio_AUC(mode,1)=trapz(SC_recon_ratio(1:mode,1))/mode;
    SC_recon_FC_r_AUC(mode,1)=trapz(SC_recon_FC_r(1:mode,1))/mode;
end
outputAUC = [(1:n_modes)',SC_recon_ratio_AUC,SC_recon_FC_r_AUC];
xlswrite('Results/CC_recon_marmoset_activity/Marmoset_SC_recon_Accuracy_AUC.xlsx',outputAUC,strcat('A2:C',num2str(n_modes+1)))

%% 1) null model SC
for i = 1: n_modes
    data(1+n_null*(i-1):n_null*i,1) = i;   
    data(1+n_null*(i-1):n_null*i,2) = SC_rewired_ratio(:,i);  
    data(1+n_null*(i-1):n_null*i,3) = SC_rewired_FC_r(:,i);  
end
xlswrite('Results/CC_recon_marmoset_activity/Marmoset_SC_recon_null_rewired.xlsx',data,strcat('A2:C',num2str(n_modes*n_null+1)))

%% 2) null model random signals
for i = 1: n_modes
    data(1+n_null*(i-1):n_null*i,1) = i;   
    data(1+n_null*(i-1):n_null*i,2) = SC_moran_ratio(:,i);  
    data(1+n_null*(i-1):n_null*i,3) = SC_moran_FC_r(:,i);  
end
xlswrite('Results/CC_recon_marmoset_activity/Marmoset_SC_recon_null_Moran.xlsx',data,strcat('A2:C',num2str(n_modes*n_null+1)))

% %=========================================================================================
% %  Calculating p-value
% % ========================================================================================
for i =1:n_modes
     % CC recon. 
    P_CC_rewired_ratio(i) = sum(CC_rewired_ratio(:,i)>CC_recon_ratio(i))/n_null;
    P_CC_rewired_FC(i) = sum(CC_rewired_FC_r(:,i)> CC_recon_FC_r(i))/n_null;
    
    P_CC_moran_ratio(i) = sum(CC_moran_ratio(:,i)>CC_recon_ratio(i))/n_null;
    P_CC_moran_FC(i) = sum(CC_moran_FC_r(:,i)> CC_recon_FC_r(i))/n_null;

    % SC recon.
    P_SC_rewired_ratio(i) = sum(SC_rewired_ratio(:,i)>SC_recon_ratio(i))/n_null;
    P_SC_rewired_FC(i) = sum(SC_rewired_FC_r(:,i)> SC_recon_FC_r(i))/ n_null;
    
    P_SC_moran_ratio(i) = sum(SC_moran_ratio(:,i)>SC_recon_ratio(i))/n_null;
    P_SC_moran_FC(i) = sum(SC_moran_FC_r(:,i)> SC_recon_FC_r(i))/ n_null;
end
P_CC_rewired_ratio(n_modes) = 0; P_CC_rewired_FC(n_modes) = 0;
P_CC_moran_ratio(n_modes) = 0; P_CC_moran_FC(n_modes) = 0;

P_SC_rewired_ratio(n_modes) = 0; P_SC_rewired_FC(n_modes) = 0;
P_SC_moran_ratio(n_modes) = 0; P_SC_moran_FC(n_modes) = 0;

null_P_CC = [P_CC_rewired_ratio;P_CC_rewired_FC;P_CC_moran_ratio;P_CC_moran_FC]';
null_P_SC = [P_SC_rewired_ratio;P_SC_rewired_FC;P_SC_moran_ratio;P_SC_moran_FC]';
xlswrite('Results/CC_recon_marmoset_activity/Marmoset_CC_SC_null_Pvalue.xlsx',null_P_CC)
xlswrite('Results/CC_recon_marmoset_activity/Marmoset_CC_SC_null_Pvalue.xlsx',null_P_SC)
