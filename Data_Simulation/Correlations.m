function CorrMat = Correlations(categorized_data,Net)
n_met = length(Net.MetabName);
n_rxn = length(Net.EnzName);

CCC = categorized_data(1:n_met,:);
FCC = categorized_data(1:n_rxn,:);

for k = 1:n_met
    for kk = 1:n_rxn  
    a = CCC(k,:);
    b = FCC(kk,:);
    %id = a ~= 0 & b ~= 0;
    %correlation(k,kk) = corr2(a(id),b(id));    
    correlation(k,kk) = corr2(a,b);
    end
end

CorrMat = correlation;
end