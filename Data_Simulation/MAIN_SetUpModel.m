%This code takes the standard flux distribution from the maranas model 
%and re-adjusts it by using metabolic flux analysis, to get closer to a
%steady state.
%If needed, metabolites and reactions can be removed. We use the cobra
%toolbox here to modify the network.
rng(1)

%set some parameters
EnsembleSize = 10000;
cutoffs = [0.0005, 2/3];


%initCobraToolbox
modelx = load('model');
modelx = modelx.model;

%%%%% setup model structure %%%%%
rxns = modelx.rxn;
mets = modelx.metab;
model.S = modelx.S;
model.mets = mets;
model.rxns = rxns;
model.Vnet = modelx.Vss;
model.Sreg = modelx.enzymeReg;

%adjust stoichiometric matrix and fluxes according to flux direction
model1 = change_directions(model);
%remove metabolites
model2 = remove_mets(model1);

%remove isolated rxns
model3 = remove_isolated_rxns(model2);

model4 = add_reg(60,model3);
[G,model] = construct_bipartite_graph(model4);

%add and remove regulation
[G,model] = add_regulation(model,G); %regulation edges hab label 2

%remove non-connected modules/ isolated modules
[G,model] = remove_nonconnected_modules(model,G);

%%%how many target genes do the regulators have
idx = G.Edges.Weight == 2;
regedges = G.Edges.EndNodes(idx,:);
mets = regedges(:,1);
reactions = regedges(:,2);
unique_mets = unique(mets);
unique_rxn = unique(reactions);
for k=1:length(unique_mets)
    count(k) = sum(unique_mets(k)==mets); 
end
for k=1:length(unique_rxn)
    count1(k) = sum(unique_rxn(k)==reactions); 
end
figure(1)
subplot(1,2,1)
histogram(count)
xlabel('#regulations per metabolite')
ylabel('#Frequency')
subplot(1,2,2)
histogram(count1);
xlabel('#metabolites that regulatate one reactions')
ylabel('#Frequency')

%% perform MCA to get flux and concentration control coefficients
Nodes = [model.mets;model.rxns];%re-define nodes
Ereg = model.Sreg;
CCC_results = cell(EnsembleSize,1);
FCC_results = cell(EnsembleSize,1);

Net.Vref = model.Vnet;
Net.EnzName = model.rxns;
Net.MetabName = model.mets;
Net.S = model.Sorig;
Net.Sreg = model.Sreg;
Net.Reversibilities = zeros(size(model.rxns));

parfor k = 1:EnsembleSize
        disp(k)
        E = setupE(model.Easo_scaffold,Ereg);
        jac = model.Sorig*diag(model.Vnet)*E';
        ew = max(real(eig(jac)));
        CCC = -1*pinv(model.Sorig*diag(model.Vnet)*E')*model.Sorig*diag(model.Vnet);
        FCC =  E'*CCC+diag(ones(length(model.Vnet),1));    
        CCC_results{k,1} = -1*CCC;
        FCC_results{k,1} = -1*FCC;
end

%%transform data / categorisation into 1,-1 etc.
[categorized_data] = categorization(CCC_results, FCC_results, cutoffs);
[minmax_data,mean_data] = minmaxing(CCC_results, FCC_results);
categorized_native = abs(categorized_data).*mean_data;


%Some Analyses
%==============================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%how many regulations can we realistically infer?
id = model.Vnet < 0.001;
fraction_low_flux = sum(id)/length(model.Vnet); %so roughly 1/3 of the fluxes have an insignificant flux
reg_low = model.Sreg(:,id);
reg_high = model.Sreg;
reg_high(:,id) = [];
total_reg = sum(model.Sreg(:) ~=0);
low_reg = sum(reg_low(:) ~=0)/total_reg;
high_reg = sum(reg_high(:) ~=0)/total_reg;
%
CCC = categorized_data(1:size(CCC_results{1},1),:);

state = cell(size(CCC,12),1);
for k = 1:size(CCC,2)
    sub_id = find(model4.S(:,k) == -1);
    sub_name = model4.mets(sub_id);
    state(k) = {max(CCC(sub_id,k))};
    if isempty(sub_name)
       state(k) = {2};
    end
end
state = cell2mat(state);
delid = state == 2;
state(delid) = [];
fraction = sum(state==1)/length(state);
fraction1 = sum(state==-1)/length(state);

%fraction of reactions where at least 1 subsrate is upregulated upon perturbation
vnet = model.Vnet;
vnet(delid) = [];
bla = state' == (vnet>0.001); % in x % of cases the substrate is upregulated when the reaction carries a high flux (above 0.001)

disp(fraction)
disp(sum(bla)/length(bla));

%find regulations which are likely to have an effect based on the
%coefficients
x = score_each_regulation(model,categorized_data,G,Nodes);
x_new = x~=0;

