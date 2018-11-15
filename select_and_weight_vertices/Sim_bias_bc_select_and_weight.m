% check if selecting and weighting the vertices you use for your analyses
% biases the laminar profiles

% TO DO
% - need to add a cross validation option to see that the u-shaped effect then
% disappears because it is due to double-dipping
% - add a proper color scale

%%
clear 
clc

% FileName = 'Sim_select_weight.tif';

nb_layers = 3;
nb_sess = 2;
nb_vertices = 2000;

layer_2_plot = 3;
session_2_plot = 1;

% Signal across layers
% Mu_cdt_1 = zeros(1,NbLayers);
% Mu_cdt_1 = -1*(1:NbLayers)/NbLayers;
% Mu_cdt_1 = [0.6970   0.6426   0.8655   1.1374   1.4821   1.8846]; % taken from empirical values
mu_cdt_1 = [0.6970   1.0014   1.8846];
mu_cdt_2 = mu_cdt_1 / 2 ;

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


%% Generate data
for iSess=1:nb_sess
    % data for vertices 1 : that respond stronger to cdt 1 than 2
    vert_1(:, :, iSess, 1) = mvnrnd(mu_cdt_1, sigma_noise, nb_vertices); %Cdt 1
    vert_1(:, :, iSess, 2) = mvnrnd(mu_cdt_2, sigma_noise, nb_vertices); %Cdt 2
    
    % data for vertices 2 : that respond stronger to cdt 2 than 1
    vert_2(:, :, iSess, 1) = mvnrnd(mu_cdt_2, sigma_noise, nb_vertices); %Cdt 1
    vert_2(:, :, iSess, 2) = mvnrnd(mu_cdt_1, sigma_noise, nb_vertices); %Cdt 1
end


%% plot responses of one group of vertices to both stimuli

close all
figure('name', ['response layer ' num2str(layer_2_plot) ' - vertices 1'])
hist(...
    cat(2, ...
    vert_1(:, layer_2_plot, session_2_plot, 1), ...
    vert_1(:, layer_2_plot, session_2_plot, 2)), 25)
xlabel('activity')
legend(...
    {'vertices prefering cdt 1: resp to cdt 1', ...
    'vertices prefering cdt 1: resp to cdt 2'})


%% Contrast conditions and plot
con_vert_1 = diff(vert_1, 1, 4); % cdt 2 - cdt 1
con_vert_2 = diff(vert_2, 1, 4);

close all
figure('name', ...
    ['contrast response layer ' num2str(layer_2_plot) ' - vertices 1 & 2'])
hist(cat(2, ...
    con_vert_1(:, layer_2_plot, session_2_plot), ...
    con_vert_2(:, layer_2_plot, session_2_plot)),  25)
xlabel('activity')
legend(...
    {'vertices prefering cdt 1: [cdt 2 - cdt 1]', ...
    'vertices prefering cdt 2: [cdt 2 - cdt 1]'})


%% Select based on average across layers
% take X top and bottom percent
X = 5;
top_perc = 1:X/100*nb_vertices*2;
bottom_perc = (nb_vertices*2-X/100*nb_vertices*2):nb_vertices*2;

% vector to know which row is for which type of vertex
vert_vect = [ones(nb_vertices, 1) ; 2*ones(nb_vertices, 1)]; 

vert = [vert_1 ; vert_2]; % concat original response
con_vert = [con_vert_1 ; con_vert_2]; % concat contrasts

mean_con_vert = mean(mean(con_vert,3), 2); % mean across layers and sessions

[B,Idx] = sort(mean_con_vert);
sort_vert_vect = vert_vect(Idx);
sort_vert = vert(Idx, :, :, :);

% only keep the top and bottom X percent
vert_vect_top_bottom = [sort_vert(top_perc,:,:,:) ; sort_vert(bottom_perc,:,:,:)];
sort_vert_vect_top_bottom = [sort_vert_vect(top_perc) ; sort_vert_vect(bottom_perc)];

% compute (prefered - not preferred) for vertices preferring condition 2
profiles_2 = diff(...
    vert_vect_top_bottom(sort_vert_vect_top_bottom==2, :, :, :), ...
    1, 4);
% compute (prefered - not preferred) for vertices preferring condition 1
profiles_1 = vert_vect_top_bottom(sort_vert_vect_top_bottom==1, :, :, 1) - ...
                vert_vect_top_bottom(sort_vert_vect_top_bottom==1, :, :, 2);


profiles =  cat(1, profiles_1, profiles_2);
mean_profiles = mean(mean(profiles,3)); % mean across sessions and vertices
std_profiles = std(mean(profiles,3)); % std across vertices of mean across sessions
            
figure('name', 'laminar profiles')
errorbar(1:nb_layers, mean_profiles, std_profiles)



return
