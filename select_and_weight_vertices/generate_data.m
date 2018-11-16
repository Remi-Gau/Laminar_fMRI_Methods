function [vert_1, vert_2] = generate_data(nb_sess, mu_cdt_1, mu_cdt_2, sigma_noise, nb_vertices)
for iSess=1:nb_sess
    % data for vertices 1 : that respond stronger to cdt 1 than 2
    vert_1(:, :, iSess, 1) = mvnrnd(mu_cdt_1, sigma_noise, nb_vertices); %#ok<*AGROW> %Cdt 1
    vert_1(:, :, iSess, 2) = mvnrnd(mu_cdt_2, sigma_noise, nb_vertices); %Cdt 2
    
    % data for vertices 2 : that respond stronger to cdt 2 than 1
    vert_2(:, :, iSess, 1) = mvnrnd(mu_cdt_2, sigma_noise, nb_vertices); %Cdt 1
    vert_2(:, :, iSess, 2) = mvnrnd(mu_cdt_1, sigma_noise, nb_vertices); %Cdt 2
end
end