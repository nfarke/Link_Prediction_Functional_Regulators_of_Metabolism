function Scores_new = score_each_regulation(model,patterns,G,Nodes)

id = G.Edges.Weight == 2;
reg_edges = G.Edges.EndNodes(id,:);

reg_edges_score = [];
%how often does the regulator change? and in which conditions?
counter = 0;
for k = 1:length(reg_edges)
    
    Reg_id = reg_edges(k,1);
    Rxn_id = reg_edges(k,2);
    Condition_id = find(patterns(Reg_id,:));
    SubstrateChanges = abs(patterns(Reg_id,Condition_id)); %number of regulator changes
    RxnChanges = abs(patterns(Rxn_id,Condition_id));
    SubRxnCoChanges(k,1) = sum(SubstrateChanges == RxnChanges);
    if isempty(SubstrateChanges)
       counter = counter+1;
    end
    %do the substrates/products of the target reaction change aswell?
    target_subs = find(model.Sorig(:,Rxn_id-length(model.mets)) == 1);
    target_prods = find(model.Sorig(:,Rxn_id-length(model.mets)) == -1);
    
    TargetSubsChanges = sum(abs(patterns(target_subs,Condition_id)),1) ~= 0;
    TargetProdChanges = sum(abs(patterns(target_prods,Condition_id)),1) ~= 0;
    Together = (TargetSubsChanges + TargetProdChanges) ~= 0;
    
    %how often does the regulator, the rxn and the substrates/products
    %change
    CoChanges = (RxnChanges + Together) == 2;
    Count = sum(CoChanges);
    reg_edges_score(k,1) = Count;
end

Scores = horzcat(SubRxnCoChanges, reg_edges_score);
Scores_new = zeros(length(id),2);
list = find(id);

for k = 1:length(list)
    pos = list(k)
    Scores_new(pos,:) = Scores(k,:);
end    
end