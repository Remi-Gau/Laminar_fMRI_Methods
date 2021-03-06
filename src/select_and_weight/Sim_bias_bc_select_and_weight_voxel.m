% check if selecting and weighting the voxels you use for your analyses
% biases the laminar profiles


%%
clear
clc
close all

nb_layers = 3;
nb_sess = 6;
nb_subj = 20;
nb_sim = 20;
nb_vertices = 750; % number of vertices in each group (preferring to respond to cdt 1 or 2)

% for plotting
plot_fig = 0;
layer_2_plot = 3;
session_2_plot = 1;

%number of top and bottom verices to take
take_top_bottom = 500;
top_perc = 1:take_top_bottom;
bottom_perc = (nb_vertices*2-(take_top_bottom-1)):nb_vertices*2;

% vector to know which row is for which type of vertex
vert_vect = [ones(nb_vertices, 1) ; 2*ones(nb_vertices, 1)];

% Signal across layers
% Mu_cdt_1 = zeros(1,NbLayers);
% Mu_cdt_1 = -1*(1:NbLayers)/NbLayers;
% Mu_cdt_1 = [0.6970   0.6426   0.8655   1.1374   1.4821   1.8846]; % taken from empirical values
mu_cdt_1 = [0.6970   1.0014   1.8846];
mu_cdt_2 = mu_cdt_1 / -1.5 ;

% Covariance matrix across layers
% A=[1 1.5 2 2 1.5 1];
% Sigma_noise = A'*A;
% Sigma_noise = [...
%     2.6485    1.9059    1.0569    0.5610    0.3431    0.3011;...
%     1.9059    2.6827    2.1034    1.1775    0.5344    0.3486;...
%     1.0569    2.1034    2.8142    2.2895    1.1996    0.5430;...
%     0.5610    1.1775    2.2895    2.9694    2.3133    1.1270;...
%     0.3431    0.5344    1.1996    2.3133    2.9294    2.1847;...
%     0.3011    0.3486    0.5430    1.1270    2.1847    3.0297]; % taken from empirical values
sigma_noise = [...
    2.6485        0.8090        0.3011;...
    0.8090        2.5907        0.8350;...
    0.3011        0.8350        3.0297]; %

for iSim = 1:nb_sim
    
    for iSubj = 1:nb_subj
        
        
        %% Generate data
        [vert_1, vert_2] = generate_data(nb_sess, mu_cdt_1, mu_cdt_2, sigma_noise, nb_vertices);
        
        
        %% plot responses of one group of vertices to both stimuli
        if plot_fig
            figure('name', ['response layer ' num2str(layer_2_plot) ' - vertices 1']) %#ok<*UNRCH>
            hist(...
                cat(2, ...
                vert_1(:, layer_2_plot, session_2_plot, 1), ...
                vert_1(:, layer_2_plot, session_2_plot, 2)), 25)
            xlabel('activity')
            legend(...
                {'vertices prefering cdt 1: resp to cdt 1', ...
                'vertices prefering cdt 1: resp to cdt 2'})
        end
        
        %% Contrast conditions and plot
        con_vert_1 = diff(vert_1, 1, 4); % cdt 2 - cdt 1
        con_vert_2 = diff(vert_2, 1, 4);
        
        ground_truth(iSim,:) = mean(mean([vert_1(:,:,:,1)-vert_1(:,:,:,2) ; con_vert_2], 3)); %#ok<*SAGROW>

        if plot_fig
            figure('name', ...
                ['contrast response layer ' num2str(layer_2_plot) ' - vertices 1 & 2'])
            hist(cat(2, ...
                con_vert_1(:, layer_2_plot, session_2_plot), ...
                con_vert_2(:, layer_2_plot, session_2_plot)),  25)
            xlabel('activity')
            legend(...
                {'vertices prefering cdt 1: [cdt 2 - cdt 1]', ...
                'vertices prefering cdt 2: [cdt 2 - cdt 1]'})
        end
        
        %% Select based on average across layers
        vert = [vert_1 ; vert_2]; % concat original response
        con_vert = [con_vert_1 ; con_vert_2]; % concat contrasts

        tval_con_vert = mean(con_vert, 3) ./ (std(con_vert, 0, 3) / size(con_vert, 3)); % t-value across sessions
        
        [sorted_tval,Idx] = sort(tval_con_vert(:));
        top_voxels = Idx(top_perc);
        bottom_voxels = Idx(bottom_perc);

        % create a matrix that says if a voxel is top (+1) or bottom (-1)
        % or neither (0)
        top_bottom_voxels = zeros(size(tval_con_vert));
        top_bottom_voxels(top_voxels) = 1;
        top_bottom_voxels(bottom_voxels) = -1;
        
        top_bottom_voxel_per_layer(iSim,:,1) = nansum(top_bottom_voxels==1);
        top_bottom_voxel_per_layer(iSim,:,2) = nansum(top_bottom_voxels==-1);
        
        % create 2 arrays one keeping only the values of the top voxels,
        % the other of the bottom voxels (everything else is nan)
        act_top_voxels = nan(size(vert,1), nb_layers, nb_sess, 2);
        act_top_voxels(repmat(top_bottom_voxels==1,[1, 1, nb_sess, 2])) = ...
            vert(repmat(top_bottom_voxels==1,[1, 1, nb_sess, 2]));
        
        act_bottom_voxels = nan(size(vert,1), nb_layers, nb_sess, 2);
        act_bottom_voxels(repmat(top_bottom_voxels==-1,[1, 1, nb_sess, 2])) = ...
            vert(repmat(top_bottom_voxels==-1,[1, 1, nb_sess, 2]));
               
        % only keep the top and bottom X percent
        act_top_bottom_voxels = [act_top_voxels ; act_bottom_voxels];
        vert_vect_top_bottom = [vert_vect ; vert_vect];
        
        mean_profiles(iSubj,:) = compute_profile_plot(...
            act_top_bottom_voxels, ...
            vert_vect_top_bottom, ...
            plot_fig);
        
        if plot_fig
            title('stim selectivity after selection: preferred - non-preferred');
        end
        
        
        %% now add weighting factor and replot
        
        weighted_act_top_bottom_voxels = act_top_bottom_voxels.* ...
            abs(repmat(tval_con_vert, [2,1, nb_sess, 2])); % weight by the unsign t-value

        weighted_mean_profiles(iSubj,:) = compute_profile_plot(...
            weighted_act_top_bottom_voxels, ...
            vert_vect_top_bottom, ...
            plot_fig);
        
        if plot_fig
            title('stim selectivity after selection and weighting: preferred - non-preferred');
        end
        
    end
    
    
    
    %% group results
    close all
    
    if plot_fig
        figure('name', 'laminar profiles', 'position', [50 50 1200 600])
        
        subplot(1,2,1)
        plot_profile(mean_profiles,0)
        title('stim selectivity after selection: preferred - non-preferred')
        
        subplot(1,2,2)
        plot_profile(weighted_mean_profiles,0)
        title('stim selectivity after selection and weighting: preferred - non-preferred')
    end
    
    [~,P_1(iSim)] = ttest(mean_profiles(:,1), mean_profiles(:,2));
    [~,P_2(iSim)] = ttest(mean_profiles(:,2), mean_profiles(:,3));
    
    
    [~,P_w_1(iSim)] = ttest(weighted_mean_profiles(:,1), weighted_mean_profiles(:,2));
    [~,P_w_2(iSim)] = ttest(weighted_mean_profiles(:,2), weighted_mean_profiles(:,3));
    
    
    sim_profiles(iSim,:,1) = mean(mean_profiles);
    sim_profiles(iSim,:,2) = mean(weighted_mean_profiles);
    
end

%% plot profiles over simulations
close all

figure('name', 'laminar profiles', 'position', [50 50 1200 600])

subplot(1,3,1)
plot_profile(ground_truth,1)
title(sprintf('stim selectivity ground truth:\n preferred - non-preferred'))

subplot(1,3,2)
plot_profile(sim_profiles(:,:,1),1)
title(sprintf('stim selectivity after selection:\n preferred - non-preferred'))

subplot(1,3,3)
plot_profile(sim_profiles(:,:,2),1)
title(sprintf('stim selectivity after selection and weighting:\n preferred - non-preferred'))

print(gcf, fullfile(pwd, 'voxel_laminar_profiles.png'), '-dpng')



figure('name', 'p curve', 'Position', [50 50 1200 600], 'Color', [1 1 1]);
subplot(2,2,1)
plot_p_curve(P_1)
title('selected: p-curve (paired t-test layer 1 & 2)')

subplot(2,2,2)
plot_p_curve(P_2)
title('selected: p-curve (paired t-test layer 2 & 3)')

subplot(2,2,3)
plot_p_curve(P_w_1)
title('selected+weighted: p-curve (paired t-test layer 1 & 2)')

subplot(2,2,4)
plot_p_curve(P_w_2)
title('selected+weighted: p-curve (paired t-test layer 2 & 3)')

print(gcf, fullfile(pwd, 'voxel_p_curves.png'), '-dpng')



figure('name', 'top - bottom voxel per layer', 'position', [50 50 1200 600])
subplot(1,2,1)
bar(mean(top_bottom_voxel_per_layer(:,:,1)))
xlabel('layers (WM ---> CSF)')
ylabel('number of top voxels')

subplot(1,2,2)
bar(mean(top_bottom_voxel_per_layer(:,:,2)))
xlabel('layers (WM ---> CSF)')
ylabel('number of bottom voxels')


print(gcf, fullfile(pwd, 'voxel_selected_voxel_per_layer.png'), '-dpng')
