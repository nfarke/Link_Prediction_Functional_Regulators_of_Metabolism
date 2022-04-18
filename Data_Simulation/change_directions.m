function model = change_directions(model)
        Vx = (model.Vnet>=0) - (model.Vnet<0);
        model.Vnet = model.Vnet.*Vx;
        model.S = model.S*diag(Vx);
end
