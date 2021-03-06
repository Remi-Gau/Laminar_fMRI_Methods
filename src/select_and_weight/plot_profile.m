function plot_profile(profiles,use_std)

if ~exist('use_std', 'var')
    use_std=1;
end
mean_profiles = nanmean(nanmean(profiles,3)); % mean across sessions and vertices
if use_std
    errorbar_profiles = nanstd(nanmean(profiles,3)); % std across vertices of mean across sessions
else
    errorbar_profiles = nanstd(nanmean(profiles,3))/(size(profiles,1)^.5);
end

MIN = nanmin(mean_profiles-errorbar_profiles)-.5;
if MIN>0
    MIN = -0.1;
end
MAX = nanmax(mean_profiles+errorbar_profiles)+.5;

hold on

errorbar(1:numel(mean_profiles), mean_profiles, errorbar_profiles)
plot([0.5 numel(mean_profiles)+.5], [0 0],  '--k')

axis([0.5 numel(mean_profiles)+.5 MIN MAX])

xlabel('layers (WM ---> CSF)')
ylabel('beta value (AU)')

end