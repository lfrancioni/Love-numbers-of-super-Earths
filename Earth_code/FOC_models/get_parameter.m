function elastic_parameter = get_parameter(r, r_in_core, r_out_core, r_mantle,...
    profile)
    if r < r_in_core
        elastic_parameter = profile(1); % inner core
    elseif r < r_out_core
        elastic_parameter = profile(2); % outer core
    elseif r < r_mantle
        elastic_parameter = profile(3); % mantle
    else
        elastic_parameter = profile(4); % crust
    end
end
