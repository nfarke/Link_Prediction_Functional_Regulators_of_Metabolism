function [G,model] = construct_bipartite_graph(model)
model.Easo_scaffold = model.S < 0; %coefficients for substrates
model.Sorig = model.S; %the original stoichiometry
model.S = model.S ~= 0;
m = size(model.S,1); %~met dimensions
n = size(model.S,2); % ~rxn dimensions

Incidence = [zeros(m,m), model.S; %incidence matrix for bipartite graph, dim: cnum + bnum x cnum + bnum,  met + gene x met + gene
        model.S', zeros(n,n)];
    
G =  graph(Incidence); %convert matrix to graph structure


%adjust Ereg
Sreg1 = model.Sreg >= -3 & model.Sreg <0; %logical negative regulation
Sreg2 = model.Sreg == -4; %logical positive regulation
model.Sreg = Sreg2 - Sreg1; %1 positive, -1 negative reg
end