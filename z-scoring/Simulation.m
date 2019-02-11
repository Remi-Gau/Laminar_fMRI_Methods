% Generate data (ignoring BOLD-response convolution) for 3 layers

b_stim_V = [0.12 0.12 0.12; ...
            0.1 0.12 0.12; ...
            0.12 0.1 0.12;...
            0.12 0.12 0.1;...
            0.10 0.12 0.10]; % param for V stim -> zscore shows maximum in middle layers


b_stim_A = 0.05;
b_H = 10;
b_const = 1;


x_V = repmat([1 1 0 0 0 1 1 0 0 0 0 0 0 0 1 1 0 1 0 1],1,1000)';
x_A = repmat([0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0],1,1000)';
x_H = repmat([1 1 1 1 1 0 1 1 1 1 0 0 0 1 1 1 0 0 0 0],1,1000)';  %% cond of no interest

X = [x_V x_A x_H ones(length(x_V),1)];

% Noise
mu = 0; sigma = 1;

Layer_scale_signal = [1 2 3]; % assume linear scaling similar to Makuerkiaga

% scaling for ERRROR ??? do we know anything about this???
Layer_scale_error = [1 1 1; ...
                    10 1 1; ...
                    1 10 1; ...
                    1 1 10; ...
                    10 1 10; ...
                    10 10 10]; 

figure; hold on

for Layer_scale_error_num = 1 : size(Layer_scale_error,1)
    for b_num = 1 : size(b_stim_V,1)
        for sim = 1 : 1

            for layer = 1 : length(Layer_scale_signal)

                % generate data for each layer
                y(:,layer) = Layer_scale_signal(layer)*(x_V* b_stim_V(b_num,layer) + x_A * b_stim_A + ones(length(x_V),1)*b_const) + x_H*b_H + Layer_scale_error(Layer_scale_error_num,layer)*normrnd(mu, sigma,length(x_V),1);

                % estimate stimulus A vs. V beta  without z-scoring
        %         beta(:,layer) = pinv(X)*y(:,layer);
                beta(:,layer) = regress(y(:,layer),X);
                stim_effect(sim,layer) =  beta(1,layer)-  beta(2,layer);

                % estimate stimulus A vs. V beta  with z-scoring
                beta_zscore(:,layer) = regress(zscore(y(:,layer)),X);
        %         beta_zscore(:,layer) = pinv(X)*zscore(y(:,layer));

                stim_effect_zscore(sim,layer) =  beta_zscore(1,layer)-  beta_zscore(2,layer);


            end


        end

        stim_mean = mean(stim_effect,1);
        stim_zscore_mean = mean(stim_effect_zscore,1);

        subplot(size(Layer_scale_error,1),size(b_stim_V,1), (Layer_scale_error_num-1)*size(b_stim_V,1)+b_num)
        plot([1 2 3], stim_mean, 'r'); hold on
        plot([1 2 3], stim_zscore_mean, 'b')
%         plot([1 2 3], b_stim_V(b_num,:), 'g')
%         plot([1 2 3], b_stim_V(b_num,:).*Layer_scale_signal, 'k')
        

    end
end
