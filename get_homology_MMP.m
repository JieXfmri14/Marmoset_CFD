function order = get_homology_MMP(raw_ROI_label,homology_ROI_label)

%% find homology index
order =zeros(max(unique(homology_ROI_label)),1);
for i = 1:length(raw_ROI_label)
    for j=1:max(unique(homology_ROI_label))
        ind = find(homology_ROI_label ==j);
        order(j) = raw_ROI_label(ind);
    end
end  