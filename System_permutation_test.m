function [net_low_sig,net_high_sig] = System_permutation_test(zX_RS,U,Vlow,Vhigh,mean_low,mean_high)

% To examine the significance of system-level concentration of low- and high-frequency signals, 
% we performed a non-parametric permutation test for each cortical class

% % ==================================================================================
% % 1) Generate CC-informed graph signal surrogates by graph spectral randomization (SR)
% % ==================================================================================
nullSurr = 1000;
data = zX_RS;
n_ROI = size(data,1);
time = size(data,2);
sub = size(data,3);
 
null_fdata=zeros(n_ROI,time,sub,nullSurr);   
for s=1:sub   
    X=data(:,:,s);
    for n=1:nullSurr
        % randomize sign of Fourier coefficients
        PHIdiag=round(rand(size(U,1),1));
        PHIdiag(PHIdiag==0)=-1; 
        PHI=diag(PHIdiag);           
        null_fdata(:,:,s,n)=U*PHI*(U'*X);   
    end
end

% % ==================================================================================
% % 2) filter low-and high-frequency signal in surrogates
% % ==================================================================================
null_X_hat_surr=zeros(n_ROI,time,sub,nullSurr);
null_X_low_surr=zeros(n_ROI,time,sub,nullSurr);
null_X_high_surr=zeros(n_ROI,time,sub,nullSurr);
null_N_low_surr=zeros(n_ROI,nullSurr,sub);
null_N_high_surr=zeros(n_ROI,nullSurr,sub);
for s=1:sub   
    for i=1:nullSurr  
        null_X_hat_surr(:,:,s,i)=U'*null_fdata(:,:,s,i);
        null_X_low_surr(:,:,s,i)=Vlow*null_X_hat_surr(:,:,s,i);     
        null_X_high_surr(:,:,s,i)=Vhigh*null_X_hat_surr(:,:,s,i);
        % norms of the filtered signals across time points
        for r=1:n_ROI
            null_N_low_surr(r,i,s)=norm(null_X_low_surr(r,:,s,i));
            null_N_high_surr(r,i,s)=norm(null_X_high_surr(r,:,s,i));
        end
    end
end
null_N_c = mean(null_N_low_surr,3);    % Average subject, roi*null
null_N_d = mean(null_N_high_surr,3);

% % ==================================================================================
% % 3) network-level
% % ==================================================================================
network_level = xlsread('marmoset_brain_template/marmoset_55Nodes_11Classes.xlsx','sheet1','D2:D56');
net_num = length(unique(network_level));
network = ["VC","PPC","PCC","LIT","AU","SS",...
     "MOT","mPFC","OFC","VLPFC","DLPFC"];

% null distribution of 11 classes
net_null_N_low = zeros(net_num,nullSurr);    
net_null_N_high = zeros(net_num,nullSurr); 
net_mean_low = zeros(net_num,1);
net_mean_high = zeros(net_num,1);
net_median_low = zeros(net_num,1);
net_median_high = zeros(net_num,1);
for i=1:net_num
        temp = find(network_level==i);
        net_null_N_low(i,:)=mean(null_N_c(temp,:),1);   % average network
        net_null_N_high(i,:)=mean(null_N_d(temp,:),1);
        net_mean_low(i)=mean(mean_low(temp));
        net_mean_high(i)=mean(mean_high(temp));  
        net_median_low(i)=median(mean_low(temp));
        net_median_high(i)=median(mean_high(temp));
end
net_null_N_low_all =zeros(nullSurr*net_num,1);
net_null_N_high_all =zeros(nullSurr*net_num,1);

% save visualization, Boxplot ordered by the median value
[~,index_cou] = sort(net_median_low,'descend');
[~,index_decou] = sort(net_median_high,'descend');

for i=1:net_num
    net_null_N_low_all(nullSurr*(i-1)+1:nullSurr*i,1)=net_null_N_low(index_cou(i),:);
    net_low_all(nullSurr*(i-1)+1:nullSurr*i,1) = network(index_cou(i));

    net_null_N_high_all(nullSurr*(i-1)+1:nullSurr*i,1)=net_null_N_high(index_decou(i),:);
    net_high_all(nullSurr*(i-1)+1:nullSurr*i,1) = network(index_decou(i));
end

path1='Result\filtered_low_high_signals';
xlswrite([path1, filesep,'network_low_compoents_all_null.xlsx'],net_null_N_low_all,'sheet1')
xlswrite([path1, filesep,'network_low_compoents_all_null.xlsx'],net_low_all,'sheet2')

xlswrite([path1, filesep,'network_high_compoents_all_null.xlsx'],net_null_N_high_all,'sheet1')
xlswrite([path1, filesep,'network_high_compoents_all_null.xlsx'],net_high_all,'sheet2')

% % ==================================================================================
% % 4)the significance of system-level concentration, one-tailed, FDR-corrected
% % ==================================================================================
%% low-frequency components
for i = 1:net_num
    if net_mean_low(i) > median(net_null_N_low(i,:))
        p_perm_low(i) = sum(net_null_N_low(i,:) > net_mean_low(i))/nullSurr;
    else
        p_perm_low(i) = 0.06;    
    end
end
FDR_low=mafdr(p_perm_low,'BHFDR', true);
net_low_sig = network(find(FDR_low<0.05));

%% high-frequency components
for i = 1:net_num
    if net_mean_high(i) > median(net_null_N_high(i,:))
        p_perm_high(i) = sum(net_null_N_high(i,:) > net_mean_high(i))/nullSurr;
    else
        p_perm_high(i) = 0.06;
    end
end
FDR_high=mafdr(p_perm_high,'BHFDR', true);
net_high_sig = network(find(FDR_high<0.05));

end