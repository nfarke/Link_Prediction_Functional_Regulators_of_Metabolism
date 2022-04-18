function [categorized] = categorization(CCC_results, FCC_results, cutoffs)
    ccc_x = cat(3, CCC_results{:});
    up = ccc_x > cutoffs(1);
    down = ccc_x < -1*cutoffs(1);
    up_fraction = sum(up,3)/size(up,3);
    down_fraction = sum(down,3)/size(down,3);
    
    up_final = up_fraction > cutoffs(2);
    down_final = down_fraction > cutoffs(2);
    categorical_ccc = up_final - down_final;
    
    fcc_x = cat(3, FCC_results{:});
    up = fcc_x > cutoffs(1);
    down = fcc_x < -1*cutoffs(1);
    up_fraction = sum(up,3)/size(up,3);
    down_fraction = sum(down,3)/size(down,3);
    
    up_final = up_fraction > cutoffs(2);
    down_final = down_fraction > cutoffs(2);
    categorical_fcc = up_final - down_final;
    
    categorized = vertcat(categorical_ccc,categorical_fcc);
end