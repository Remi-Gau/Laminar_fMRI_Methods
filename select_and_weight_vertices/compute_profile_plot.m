function mean_profiles = compute_profile_plot(vect_act, vect_id, plot_fig)

% compute (prefered - not preferred) for vertices preferring condition 2
profiles_2 = diff(...
    vect_act(vect_id==2, :, :, :), ...
    1, 4);
% compute (prefered - not preferred) for vertices preferring condition 1
profiles_1 = vect_act(vect_id==1, :, :, 1) - ...
    vect_act(vect_id==1, :, :, 2);

profiles =  cat(1, profiles_1, profiles_2);

mean_profiles = mean(mean(profiles,3)); % mean across sessions and vertices

if plot_fig
    figure('name', 'laminar profile')
    plot_profile(profiles);
end

end