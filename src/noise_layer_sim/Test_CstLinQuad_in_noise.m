% script that generates null data to check if the cst, lin and quad
% component of the laminar GLM are on average equal to 0 and give a flat
% p-curve

clear; clc; close all

% Start_dir = fullfile('D:\Dropbox','PhD','Experiments','Laminar_fMRI_Methods');
% addpath(genpath(fullfile(Start_dir,'code')))
% Get_dependencies('D:\Dropbox/')

IID = 1; % in case we want the data across layers to be iid or not
Do_perm = 0; % perform t-test or sign permutation test
Do_Var_Expl = 0; % compute percent variance explained

NbLayers = 6;
NbVertices = 1000; %5000
NbSubj = 10;
NbSess = 20;

NbSim = 1000; %10000
PrintOutEvery = 10^3;

Mu = zeros(1,NbLayers); % mean value

if ~IID
    NoiseSuffix = ' - layers not iid';
    % covariance matrix to generate data: taken from real data
    Sigma_noise = [...
        2.6485    1.9059    1.0569    0.5610    0.3431    0.3011;...
        1.9059    2.6827    2.1034    1.1775    0.5344    0.3486;...
        1.0569    2.1034    2.8142    2.2895    1.1996    0.5430;...
        0.5610    1.1775    2.2895    2.9694    2.3133    1.1270;...
        0.3431    0.5344    1.1996    2.3133    2.9294    2.1847;...
        0.3011    0.3486    0.5430    1.1270    2.1847    3.0297];
else
    Sigma_noise = eye(NbLayers);
    NoiseSuffix = ' - layers iid';
end

% Design matrix to use for the laminar GLM
DesMat = (1:NbLayers)-mean(1:NbLayers);
DesMat = [ones(NbLayers,1) DesMat' (DesMat.^2)'];
DesMat = spm_orth(DesMat);
X=repmat(DesMat,NbSess,1);

if Do_perm
    % Cartesian product to list all the possible sign permutation with that
    % many subjects.
    for iSub=1:NbSubj
        sets{iSub} = [-1 1];
    end
    [a, b, c, d, e, f, g, h, i, j] = ndgrid(sets{:});
    ToPermute = [a(:), b(:), c(:), d(:), e(:), f(:), g(:), h(:), i(:), j(:)];
    clear a b c d e f g h i j
end

%%
t = tic;
T = [];
for iSim = 1:NbSim
    
    % Just to diplay how time we have left to run all the simulations
    if mod(iSim,PrintOutEvery)==0
        
        T(end+1)=toc(t);
        
        sec = round((NbSim-iSim)*mean(T/PrintOutEvery));
        hrs = floor(sec/3600);
        min = floor(mod(sec,3600)/60);
        
        fprintf('Simulation %i\n', iSim)
        fprintf(1,'Avg time elapsed / simulation = %0.3f secs ; ETA = %i hrs %i min\n',...
            mean(T/PrintOutEvery), hrs, min);
        
        t = tic;
    end
    
    % Fits the GLM for each subject
    for iSubj = 1:NbSubj
        
        % Generate data for each session
        for iSess=1:NbSess
            if ~IID
                Dist(:,:,iSess) = mvnrnd(Mu, Sigma_noise, NbVertices); %#ok<*SAGROW>
            else
                Dist(:,:,iSess) = randn(NbVertices,NbLayers); %#ok<*SAGROW>
            end
        end
        
        % Change or adapt dimensions for GLM
        Y = shiftdim(Dist,1);
        Y = reshape(Y, [size(Y,1)*size(Y,2), size(Y,3)] );
        
        % GLM for each vertex/session and take the mean of that over
        % session
        B = pinv(X)*Y;
        betas(iSubj,1:size(DesMat,2)) = mean(B,2);
        
        if Do_Var_Expl
            % Variance explained under the full model
            Y_hat = X * B;
            Res  = Y - Y_hat;
            
            SS_total = norm(Y - mean(Y)).^2;
            SS_effect = norm(Y_hat - mean(Y_hat)).^2;
            R2 = SS_effect / SS_total;
            clear Y_hat Res SS_effect
            
            % Variance explained by each regressor of the matrix
            % (http://www.sbirc.ed.ac.uk/cyril/glm/GLM_lectures.html#5)
            % explained variance --> squared semi partial correlation
            for iReg = 1:size(DesMat,2)
                X_red = X;
                X_red(:,iReg) = []; % reduced model all minus 1 regressor
                
                B_red = pinv(X_red)*Y;
                
                Y_hat_red = X_red*B_red;
                
                SS_effect_red = norm(Y_hat_red-mean(Y_hat_red)).^2;
                R2_reduced = SS_effect_red / SS_total;
                
                Semi_Partial_corr_coef(iReg,iSubj) = R2 - R2_reduced;
                
                clear Y_hat_red B_red Res SS_effect_red R2_reduced
            end
            
        end
        
    end
    
    % P value using t-test
    [h,p]=ttest(betas, 0, 'tail', 'both');
    Results_ttest(iSim,:) = p;
    
    if Do_perm
        
        % Do sign permutation tests accross subjects
        for ibeta = 1:size(betas,2)
            tmp = repmat(betas(:,ibeta)',size(ToPermute,1),1);
            Perms = mean(ToPermute.*tmp,2);
            p(ibeta) = sum( abs(Perms) > abs(mean(betas(:,ibeta))) ) / numel(Perms) ;
            Nulls(:,ibeta,iSim)=Perms;
        end
        
        Results_perm(iSim,:) = p;
        
    end
    
end

save(fullfile(Start_dir,'results','p-curve_profiles',...
    ['simulation_profiles_cst_lin_quad' NoiseSuffix ' - ' datestr(now, 'yyyy_mm_dd_HH_MM')   '.mat']))


%% Plot distribution across subjects of the variance explained by each regressor
if Do_Var_Expl
    
    subplot(311)
    hist(Semi_Partial_corr_coef(1,:),50)
    title('R_{full}^2 - R_{reduced_{Cst}}^2')
    
    subplot(312)
    hist(Semi_Partial_corr_coef(2,:),50)
    title('R_{full}^2 - R_{reduced_{Lin}}^2')
    
    subplot(313)
    hist(Semi_Partial_corr_coef(3,:),50)
    title('R_{full}^2 - R_{reduced_{Quad}}^2')
    
    
    p=mtit(['Squared semi partial corr coeff' NoiseSuffix],...
        'fontsize',14,...
        'xoff',0,'yoff',.05);
    
end

%% plot p curves
NbRow = 3;
NbCol = 1;
Subplot = 1:3;
if Do_perm
    NbCol = 2;
    Subplot = [1 3 5 2 4 6];
end

if NbSim>1

    NbBins = 100;
    
    close all
    
    figure('name', ['Simulation Cst Lin Quad' NoiseSuffix], 'Position', [100, 100, 1000, 700], 'Color', [1 1 1]);
    
    subplot(NbRow,NbCol,Subplot(1))
    hist(Results_ttest(:,1),NbBins);
    H(1,1:NbBins) = hist(Results_ttest(:,1),NbBins);
    ylabel('Constant')
    title('T-test')
    
    subplot(NbRow,NbCol,Subplot(2))
    hist(Results_ttest(:,2),NbBins);
    H(2,1:NbBins) = hist(Results_ttest(:,2),NbBins);
    ylabel('Linear')
    
    subplot(NbRow,NbCol,Subplot(3))
    hist(Results_ttest(:,3),NbBins);
    H(3,1:NbBins) = hist(Results_ttest(:,3),NbBins);
    ylabel('Quadratic')
    
    if Do_perm
        subplot(NbRow,NbCol,Subplot(4))
        hist(Results_perm(:,1),NbBins);
        H(4,1:NbBins) = hist(Results_perm(:,1),NbBins);
        title('Permutation test')
        
        subplot(NbRow,NbCol,Subplot(5))
        hist(Results_perm(:,2),NbBins);
        H(5,1:NbBins) = hist(Results_perm(:,2),NbBins);
        
        subplot(NbRow,NbCol,Subplot(6))
        hist(Results_perm(:,1),NbBins);
        H(6,1:NbBins) = hist(Results_perm(:,3),NbBins);
    end
    
    for i=1:numel(Subplot)
        subplot(NbRow,NbCol,Subplot(i))
        hold on
        
        plot([.05 .05], [0 max(H(:))], 'r', 'linewidth', 2)
        
        [x,y] = ind2sub([3,2],i);
        switch y
            case 1
                tmp=Results_ttest;
            case 2
                if Do_perm
                    tmp=Results_perm;
                end
        end
        text(.7, max(H(:)), sprintf('p(p<.05)=%.04f',mean(tmp(:,x)<.05)))
        
        set(gca, 'xtick', 0:.1:1,'xticklabel', 0:.1:1)
        xlabel('p-value')
        axis([0 1 0 max(H(:))*1.1])
    end
    
end
