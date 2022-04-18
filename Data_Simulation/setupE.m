function [E] = setupE(Easo_scaffold,Ereg)
    LB = log10(0.001); %define lower bound
    UB = log10(1); %define upper bound
    id_easo = find(Easo_scaffold);
    Easo1 = zeros(size(Easo_scaffold));
    Easo1(id_easo) = 10.^(LB+(UB-LB)*rand(1,length(id_easo)));    
    Easo = Easo1;
    LB = log10(1); %define lower bound
    UB = log10(4); %define upper bound
    Ereg(find(Ereg)) = Ereg(find(Ereg)).*10.^(LB+(UB-LB)*rand(1,length(find(Ereg))))';
    
    E = Easo + Ereg;
end