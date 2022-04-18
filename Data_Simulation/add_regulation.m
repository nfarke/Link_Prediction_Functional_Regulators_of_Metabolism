function [G,model] = add_regulation(model,G)

%add regulation
[met_id,rxn_id] = ind2sub(size(model.Sreg),find(model.Sreg));
mode = model.Sreg(find(model.Sreg));
reg_edges = [met_id,length(model.mets)+rxn_id,mode];

%remove substrate/product regulation
[~,id1,id2] = intersect(G.Edges.EndNodes,reg_edges(:,1:2),'rows');
reg_edges(id2,:) = [];

Ereg_new = zeros(size(model.Sreg));
linearInd = sub2ind(size(Ereg_new), reg_edges(:,1), reg_edges(:,2)-length(model.mets));
Ereg_new(linearInd) = reg_edges(:,3);
model.Sreg = Ereg_new;
G = addedge(G,reg_edges(:,1),reg_edges(:,2),2*ones(length(reg_edges),1));
end