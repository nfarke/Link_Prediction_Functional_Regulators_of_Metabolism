function model = add_reg(num_reg, model)
Sreg = model.Sreg;
%find metabolites that are not yet involved in regulation
for k = 1:size(Sreg,1)
    if isempty(find(Sreg(k,:)))
        val = 1;
    else
        val = 0;
    end
    met_ids(k) = logical(val);
end
met_ids = find(met_ids);

%find reactions that are not yet involved in regulation
for k = 1:size(Sreg,2)
    if isempty(find(Sreg(:,k)))
        val = 1;
    else
        val = 0;
    end
    rxn_ids(k) = logical(val);
end
rxn_ids = find(rxn_ids);
significant_flux_ids = model.Vnet(rxn_ids) > 0.001;
rxn_ids = rxn_ids(significant_flux_ids);

%randomly choose regulator and reactions
metx = randsample(met_ids,num_reg,0);
rxnx = randsample(rxn_ids,num_reg,0);

for j = 1:num_reg
    mety = metx(j);
    rxny = rxnx(j);
    model.Sreg(mety,rxny) = -1; %negative regulation
    if j > 0.9*num_reg
    model.Sreg(mety,rxny) = -4;    
    end
end

end