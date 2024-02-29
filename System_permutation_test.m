function [class_low_sig,class_high_sig] = System_permutation_test(zX_RS,U,Vlow,Vhigh,mean_low,mean_high)

% To examine the significance of system-level concentration of low- and high-frequency signals, 
% we performed a non-parametric permutation test for each cortical class

% % ==================================================================================
% % 1) Generate CC-informed graph signal surrogates by graph spectral randomization (SR)
% % ==================================================================================
n_null = 1000;
data = zX_RS;
n_ROI = size(data,1);
time = size(data,2);
sub = size(data,3);
 
null_fdata=zeros(n_ROI,time,sub,n_null);   
for s=1:sub   
    X=data(:,:,s);
    for n=1:n_null
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
null_X_hat_surr=zeros(n_ROI,time,sub,n_null);
null_X_low_surr=zeros(n_ROI,time,sub,n_null);
null_X_high_surr=zeros(n_ROI,time,sub,n_null);
null_N_low_surr=zeros(n_ROI,n_null,sub);
null_N_high_surr=zeros(n_ROI,n_null,sub);
for s=1:sub   
    for i=1:n_null  
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
% % 3) class-level
% % ==================================================================================
class_level = xlsread('marmoset_brain_template/marmoset_55Nodes_11Classes.xlsx','Sheet1','D2:D56');
n_net = length(unique(class_level));
network = ["VC","PPC","PCC","LIT","AU","SS",...
     "MOT","mPFC","OFC","VLPFC","DLPFC"];

% null distribution of 11 classes
net_null_low = zeros(n_net,n_null);    
net_null_high = zeros(n_net,n_null); 
net_mean_low = zeros(n_net,1);
net_mean_high = zeros(n_net,1);
net_median_low = zeros(n_net,1);
net_median_high = zeros(n_net,1);
for i=1:n_net
        temp = find(class_level==i);
        net_null_low(i,:)=mean(null_N_c(temp,:),1);   % average network
        net_null_high(i,:)=mean(null_N_d(temp,:),1);
        net_mean_low(i)=mean(mean_low(temp));
        net_mean_high(i)=mean(mean_high(temp));  
        net_median_low(i)=median(mean_low(temp));
        net_median_high(i)=median(mean_high(temp));
end
net_null_low_all =zeros(n_null*n_net,1);
net_null_high_all =zeros(n_null*n_net,1);

% save visualization, Boxplot ordered by the median value
[~,index_cou] = sort(net_median_low,'descend');
[~,index_decou] = sort(net_median_high,'descend');

for i=1:n_net
    net_null_low_all(n_null*(i-1)+1:n_null*i,1)=net_null_low(index_cou(i),:);
    net_low_all(n_null*(i-1)+1:n_null*i,1) = network(index_cou(i));

    net_null_high_all(n_null*(i-1)+1:n_null*i,1)=net_null_high(index_decou(i),:);
    net_high_all(n_null*(i-1)+1:n_null*i,1) = network(index_decou(i));
end

xlswrite('Results\filtered_low_high_signals\class_low_compoents_null.xlsx',net_null_low_all,'Sheet1')
xlswrite('Results\filtered_low_high_signals\class_low_compoents_null.xlsx',net_low_all,'Sheet2')

xlswrite('Results\filtered_low_high_signals\class_high_compoents_null.xlsx',net_null_high_all,'Sheet1')
xlswrite('Results\filtered_low_high_signals\class_high_compoents_null.xlsx',net_high_all,'Sheet2')

% % ==================================================================================
% % 4)the significance of system-level concentration, one-tailed, FDR-corrected
% % ==================================================================================
%% low-frequency components
for i = 1:n_net
    if net_mean_low(i) > median(net_null_low(i,:))
        p_perm_low(i) = sum(net_null_low(i,:) > net_mean_low(i))/n_null;
    else
        p_perm_low(i) = 0.06;    
    end
end
FDR_low=mafdr(p_perm_low,'BHFDR', true);
class_low_sig = network(find(FDR_low<0.05));

%% high-frequency components
for i = 1:n_net
    if net_mean_high(i) > median(net_null_high(i,:))
        p_perm_high(i) = sum(net_null_high(i,:) > net_mean_high(i))/n_null;
    else
        p_perm_high(i) = 0.06;
    end
end
FDR_high=mafdr(p_perm_high,'BHFDR', true);
class_high_sig = network(find(FDR_high<0.05));

end