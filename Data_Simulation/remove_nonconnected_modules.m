function [G,model] = remove_nonconnected_modules(model,G)
distx = distances(G);
inf_ids = find(isinf(distx(:,1))); 
G = rmnode(G,inf_ids);
Nodes = [model.mets;model.rxns];
Nodes2Remove = Nodes(inf_ids);

%remove them also in the model
met2remove = [];
rxn2remove = [];
for k = 1:length(Nodes2Remove)
    met_val = find(strcmp(Nodes2Remove{k},model.mets));
    if ~isempty(met_val)
        met2remove(k,1) = met_val;
    end
    rxn_val = find(strcmp(Nodes2Remove{k},model.rxns));
    
    if ~isempty(rxn_val)
        rxn2remove(k,1) = rxn_val;
    end
end

model = removeMetabolites(model,model.mets(met2remove),false);
model = removeRxns(model,model.rxns(rxn2remove(rxn2remove~=0)),false);
end