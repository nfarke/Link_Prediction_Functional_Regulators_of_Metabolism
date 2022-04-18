function model = remove_isolated_rxns(model)
for k = 1:size(model.S,2)
    if isempty(find(model.S(:,k)))
        idr(k) = 1;
    else
        idr(k) = 0;
    end
end
model = removeRxns(model,model.rxns(logical(idr)),false);

end