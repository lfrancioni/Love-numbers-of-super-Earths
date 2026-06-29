%-------------------------------------------------------------------------%
% --- Extract the value of an elastic parameter at radial distance -------%
% --- r inside the differentiated planet ---------------------------------%
%-------------------------------------------------------------------------%

function elastic_parameter = get_parameter(r, r_core, r_mantle,...
    profile)
%at the passage from a layer to another, the elastic parameter is the one
%of the outer layer (due to < sign)
    if r < r_core
        elastic_parameter = profile(1); % core
    elseif r < r_mantle
        elastic_parameter = profile(2); % mantle
    else
        elastic_parameter = profile(3); % crust
    end
end
