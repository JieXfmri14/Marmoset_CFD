function P = verify_stron_con(W)
% Verify the directed graph with strong connected
n_ROI = size(W,2);  

P = zeros(n_ROI,n_ROI);
for k=1:n_ROI
    C1 = W^k;
    P = P + C1;
end
P = (P~=0);

if P ~=1
    sprintf('This is a not strongly-connected graph')
end

end

