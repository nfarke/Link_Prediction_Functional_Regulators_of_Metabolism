# Link_Prediction_Functional_Regulators_of_Metabolism

The goal of this project was to set up a link prediction algorithm to predict allosteric regulation from metabolite data. There are many approaches to infer regulations
however these approaches are either based on mechanistic models or purely data driven disregarding available structural information. Here, we wanted to close this gap by
setting up a link prediction algorithm that includes metabolite and flux data in conjunction with structural network information. Since there is no experimental dataset for this purpose
we used a mathematical model of E.coli to simulate our data.

Simulating the data:
- Wie used the published Maranas model of E.coli
- We defined a steady state with metabolic flux analysis
- We set up the model as bipartite graph where metabolites and reactions are nodes. Substrate-Product relationships are the backbone of the graph. On top, regulations are added
  and labeled
- Because the number of known allosteric regulations are scarce and often single metabolites affect multiple reactions and vice-versa, we added 60 artificial regulations 
  where source and target were unique (later refered to unique regulations)
- We then used Metabolic Control Analysis to simulate Concentration Control Coefficients and Flux Control Coefficiens. These coeffients contain information about changes in metabolites
  or fluxes. Usually, the substrates of perturbed reactions accumulated.

Link Prediction in PyTorch:
- We chose PyTorch for its flexibility and at the time of the project PyTorch Geometric and geometric deep learning was just emerging
- The machine learning task was a binary classification: Either an edge is a regulator or not
- Since we were lacking negative examples we used the principle of negative sampling to label a certain amount of edges which carry a negative label
- We then used GCNConv with a single graph convolutional step. Here, each node in the graph learns about its neighborhood and by repeating this step for each node information
  is propagated
- We also tested how the link prediction performs when no structural information is provided (GCNConv is replaced by a simple linear layer)
- We then performed feature engineering and tested what kind of input information yields the best results (e.g. only metabolite data, metabolite + flux data, min-max scaling,
  categorized data). In this code we used categorized binary data and "No" data (only network information)
- We trained our model over several epochs using edge dropout as regularization technique
- cross-validation was used on randomized sets of edges
- We then tested several validation metrics: Because the dataset is imbalanced metrics like balanced accuracy, Precision and Recall are important. It is important to avoid false negatives
  as the number of negative links outweights the positive cases. Accuracy itself is meaningless on it´s own because a high accuracy can be obtained by classifying all cases as negatives for example.
-------------------
I designed this project. Martin Lempp then joined me and we setup the codes together.
