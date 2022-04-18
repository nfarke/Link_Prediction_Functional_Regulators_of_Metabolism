function [minmax_data, mean_data] = minmaxing(CCC_results,FCC_results)
CCC_mean = mean(cat(3, CCC_results{:}), 3);
FCC_mean = mean(cat(3, FCC_results{:}), 3);

for k=1:size(CCC_mean,1)
    CCC_new(k,:) = rescale(CCC_mean(k,:),-1,1);
end

for k=1:size(FCC_mean,1)
    FCC_new(k,:) = rescale(FCC_mean(k,:),-1,1);
end
minmax_data = [CCC_new;FCC_new];

mean_data = [CCC_mean;FCC_mean];
end