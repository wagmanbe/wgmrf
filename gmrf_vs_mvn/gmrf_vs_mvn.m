% BMW
% May 18 2020
clear all; close all
addpath('export_fig-master')
warning('off', 'MATLAB:MKDIR:DirectoryExists');
%%%% Parameters for the user to set
% Params for the obs cov. matrix. "Sigma" 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xfield_c = 0.5;                    % [0 to 0.9] cross-field correlation at the grid point
infield_sc = 0.9;                   % [0 to 1] scale the spatial covariance within fields.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note, one can tell that when using very small values for these
% parameters, GMRF creates an artifical pattern based on Q that imposes
% nonexistent relationships. 


% Load pre-calculated Q operator for a 3x3 grid, and pre-calculated alpha
% for that grid. 
Q      = dlmread('Q/Q_3x3.txt'); % works in 2016b
alpha  = dlmread('Q/alpha_3x3.txt');

%%% Part 1: Generate an invertible sigma for 3x3 grid, 2 fields, from which we will draw vectors of obs. %%%%%
[x1 x2] = deal( 0:1:2 );
[X1,X2] = meshgrid(x1,x2);
[sigma_mvn_obs, inv_sigma_mvn_obs] = make_sigma(x1, x2, xfield_c, infield_sc);

%%% Part 2: Draw 2 fields of obs from Sigma %%%%%
obs_mu = zeros(1,18); % the means
[obs1, obs2]= draw_from_mvn( x1, x2, obs_mu, sigma_mvn_obs, 100);
obs1_climo = mean( obs1, 3); obs2_climo = mean( obs2, 3);

% Calculate the spatial average covariance and its inverse.
% This is as simple as computing the covariance between vectors 
% This calculation is a spatial average because obs1(:) and obs2(:) are vectors comprised
% of timeseries from all the grid points appended together. It is the
% same thing as taking the covariance at each point and then averaging
% that in space, which is what Alvaro N-S and BMW do in the real code.
S_gmrf = cov( obs1(:), obs2(:)); % 

%%% Part 3: Generate model covariance matrix %%%%%
% Can be the same or different as the obs covariance matrix.
[sigma_mvn_mod, inv_sigma_mvn_mod] = make_sigma(x1, x2, xfield_c, infield_sc);





%%% Part 4: Draw 2 fields of model output from model covariance matrix %%%%%
n = 100;
% Test 1: Unbiased means for calculating effective degrees of freedom (ke) and scalar P_q
figdir = 'dof';
[mod1, mod2]= draw_from_mvn( x1, x2, obs_mu, sigma_mvn_mod, n);
plot_random_grid_point_corr( x1,x2, mod1,mod2)
error1 = mod1-obs1_climo; error2 = mod2-obs2_climo;

% Compute costs with p_q initialized as 1, because we do not yet know the degrees of freedom.  
[E_mvn, E_gmrf, E_trad] = compute_cost( error1, error2, Q, alpha,sigma_mvn_obs, S_gmrf, [1 1 1], figdir);

% Compute effective degrees of freedom ke.
% Note, this is the dof for the entire argument to the exponential:
% 1/2 * [(m-d)' sigma^-1 (m-d)], which is distributed A * 1/2 * Chi^2(ke) (A is an unknown constant)
% The 1/2 term comes from the 1/2 out in front of the [(m-d) sigma^-1 (m-d) ]
ke_mvn =  (2 * mean(E_mvn )^2 ) / var(E_mvn)  % Eqn. 19, Jackson and Huerta 2016 
ke_gmrf = (2 * mean(E_gmrf)^2 ) / var(E_gmrf)
ke_trad = (2 * mean(E_trad)^2 ) / var(E_trad)

pq_mvn =  sqrt(0.5 * ke_mvn) / std(E_mvn);         % Eqn. 20, Jackson and Huerta 2016
pq_gmrf = sqrt(0.5 * ke_gmrf) / std(E_gmrf);        
pq_trad = sqrt(0.5 * ke_trad) / std(E_trad); 
p_q = [pq_mvn pq_gmrf pq_trad]

% Test 0. Perfect models, normalized by p_q. 
figdir = 'perfect'
[mod1, mod2]= draw_from_mvn( x1, x2, obs_mu, sigma_mvn_mod, n);
error1 = mod1-obs1_climo;
error2 = mod2-obs2_climo;
% provide real p_q to scale errors
[cost_mvn_mean_perfect, cost_gmrf_mean_perfect, cost_trad_mean_perfect]...
    = compute_cost( error1, error2, Q, alpha,sigma_mvn_obs, S_gmrf, p_q, figdir);

% Test 1. Bias in means in both fields. 
% Shift the mean by two standard deviations. 
figdir = 'uniform_bias'
std_devs = sqrt(diag( sigma_mvn_mod ))';
[mod1, mod2]= draw_from_mvn( x1, x2, obs_mu + 2*std_devs, sigma_mvn_mod, n);
error1 = mod1-obs1_climo;
error2 = mod2-obs2_climo;
% provide real p_q to scale errors
[cost_mvn_mean_bias, cost_gmrf_mean_bias, cost_trad_mean_bias]...
    = compute_cost( error1, error2, Q, alpha,sigma_mvn_obs, S_gmrf, p_q, figdir);


% Test 2. Bias in means one field but not the other. 
figdir = 'bias_one_field'
std_devs = sqrt(diag( sigma_mvn_mod ))';
[mod1, mod2]= draw_from_mvn( x1, x2, obs_mu + [zeros(1,9) ones(1,9)].* (std_devs * 2) , sigma_mvn_mod, n);
error1 = mod1-obs1_climo;
error2 = mod2-obs2_climo;
[cost_mvn_mean_bias_one_field, cost_gmrf_mean_bias_one_field, cost_trad_mean_bias_one_field] ...
    = compute_cost( error1, error2, Q, alpha, sigma_mvn_obs, S_gmrf, p_q, figdir);


% Test 3. Bias in the covariance between fields 
figdir = 'bias_in_cross_field_covariance'
[sigma_mvn_mod_bias_c, inv_sigma_mvn_mod_bias_c] = make_sigma(x1, x2, 0.25*xfield_c, infield_sc);
[mod1, mod2]= draw_from_mvn( x1, x2, obs_mu , sigma_mvn_mod_bias_c, n);
error1 = mod1-obs1_climo;
error2 = mod2-obs2_climo;
[cost_mvn_bias_xfield_cov, cost_gmrf_bias_xfield_cov, cost_trad_bias_xfield_cov]...
    = compute_cost( error1, error2, Q, alpha, sigma_mvn_obs, S_gmrf, p_q, figdir);

% Test 4. Bias in the covariance within fields
figdir = 'bias_magnitude_of_spatial_covariance'
[sigma_mvn_mod, inv_sigma_mvn_mod] = make_sigma(x1, x2, xfield_c, 0.25 * infield_sc);
[mod1, mod2]= draw_from_mvn( x1, x2, obs_mu , sigma_mvn_mod, n);
error1 = mod1-obs1_climo;
error2 = mod2-obs2_climo;
[cost_mvn_bias_infield_cov, cost_gmrf_bias_infield_cov, cost_trad_bias_infield_cov]...
    = compute_cost( error1, error2, Q, alpha, sigma_mvn_obs, S_gmrf, p_q, figdir);

% Test 5. Bias in the covariance within and accross fields
figdir = 'bias_magnitude_of_spatial_and_xfield_covariance'
[sigma_mvn_mod, inv_sigma_mvn_mod] = make_sigma(x1, x2, 0.25 * xfield_c, 0.25 * infield_sc);
[mod1, mod2]= draw_from_mvn( x1, x2, obs_mu , sigma_mvn_mod, n);
error1 = mod1-obs1_climo;
error2 = mod2-obs2_climo;
[cost_mvn_bias_inout_cov, cost_gmrf_bias_inout_cov, cost_trad_bias_inout_cov]...
    = compute_cost( error1, error2, Q, alpha, sigma_mvn_obs, S_gmrf, p_q, figdir);

% Test 6. Bias in the mean and in the covariance within and accross fields
figdir = 'bias_magnitude_of_mean_and_spatial_and_xfield_covariance'
[sigma_mvn_mod, inv_sigma_mvn_mod] = make_sigma(x1, x2, 0.25 * xfield_c, 0.25 * infield_sc);
std_devs = sqrt(diag( sigma_mvn_mod ))';
[mod1, mod2]= draw_from_mvn( x1, x2, obs_mu + 2 * std_devs , sigma_mvn_mod, n);
error1 = mod1-obs1_climo;
error2 = mod2-obs2_climo;
[cost_mvn_bias_all, cost_gmrf_bias_all, cost_trad_bias_all]...
    = compute_cost( error1, error2, Q, alpha, sigma_mvn_obs, S_gmrf, p_q, figdir);

% Test 7. Apply a spatial gradient pattern of bias to one field.
% The bias will be 1 std in the first column, 0.5 std in the second column.
figdir = 'bias_gradient_field'
[sigma_mvn_mod, inv_sigma_mvn_mod] =make_sigma(x1, x2, xfield_c, infield_sc);
[mod1, mod2]= draw_from_mvn( x1, x2, obs_mu, sigma_mvn_mod, n);
stddevs=std(mod1,1,3);
for i=1:length(mod1)
    orig=mod1(:,:,i);
    bias = stddevs .* repmat( [1 0.5 0], [3,1]); % bias is largest in the left column
    mod1(:,:,i) = orig + abs(bias);
end
error1 = mod1-obs1_climo;
error2 = mod2-obs2_climo;
[cost_mvn_bias_grad_one_field, cost_gmrf_bias_grad_one_field, cost_trad_bias_grad_one_field] ...
    = compute_cost( error1, error2, Q, alpha, sigma_mvn_obs, S_gmrf, p_q, figdir);



% Test 8. Apply a consistent spatial gradient pattern of bias to both fields
figdir = 'bias_gradient_two_field'
[sigma_mvn_mod, inv_sigma_mvn_mod] = make_sigma(x1, x2, xfield_c, infield_sc);
[mod1, mod2]= draw_from_mvn( x1, x2, obs_mu, sigma_mvn_mod, n);
stddevs1=std(mod1,1,3);
stddevs2=std(mod2,1,3);
for i=1:length(mod1)
    orig1=mod1(:,:,i); orig2 = mod2(:,:,i); 
    bias1 = stddevs1 .* repmat( [1 0.5 0], [3,1]); % bias is largest in the left column
    bias2 = stddevs2 .* repmat( [1 0.5 0], [3,1]); % bias is largest in the left column
    mod1(:,:,i) = orig1 + abs(bias1);
    mod2(:,:,i) = orig2 + abs(bias2);

end
error1 = mod1-obs1_climo;
error2 = mod2-obs2_climo;
[cost_mvn_bias_grad_both_field, cost_gmrf_bias_grad_both_field, cost_trad_bias_both_field] ...
    = compute_cost( error1, error2, Q, alpha, sigma_mvn_obs, S_gmrf, p_q, figdir);

% Test 9. Apply different spatial gradient patterns of bias to both fields
figdir = 'bias_gradient_two_field_conflicting'
[sigma_mvn_mod, inv_sigma_mvn_mod] = make_sigma(x1, x2, xfield_c, infield_sc);
[mod1, mod2]= draw_from_mvn( x1, x2, obs_mu, sigma_mvn_mod, n);
stddevs1=std(mod1,1,3);
stddevs2=std(mod2,1,3);
for i=1:length(mod1)
    orig1=mod1(:,:,i); orig2 = mod2(:,:,i); 
    bias1 = stddevs1 .* repmat( [1 0.5 0], [3,1]); % bias is largest in the first column
    bias2 = stddevs2 .* repmat( [0 0.5 1]', [1,3]); % bias is largest in the last row
    mod1(:,:,i) = orig1 + abs(bias1);
    mod2(:,:,i) = orig2 + abs(bias2);

end
error1 = mod1-obs1_climo;
error2 = mod2-obs2_climo;
[cost_mvn_bias_grad_both_field_conf, cost_gmrf_bias_grad_both_field_conf, cost_trad_bias_grad_both_field_conf] ...
    = compute_cost( error1, error2, Q, alpha, sigma_mvn_obs, S_gmrf, p_q, figdir);


% Test 10. Apply a uniform 2 sigma bias to both fields, compute the cost before and after.
figdir = 'uniform_before_and_after'
[sigma_mvn_mod, inv_sigma_mvn_mod] = make_sigma(x1, x2, xfield_c, infield_sc);
[mod1, mod2]= draw_from_mvn( x1, x2, obs_mu, sigma_mvn_mod, n);
stddevs1=std(mod1,1,3);
stddevs2=std(mod2,1,3);
for i=1:length(mod1)
    orig1=mod1(:,:,i); orig2 = mod2(:,:,i); 
    bias1 = 2 * stddevs1 .* repmat( [1 1 1], [3,1]); 
    bias2 = 2 * stddevs2 .* repmat( [1 1 1], [3,1]); 
    mod1_biased(:,:,i) = orig1 + abs(bias1);
    mod2_biased(:,:,i) = orig2 + abs(bias2);  
end
error1 = mod1-obs1_climo;
error2 = mod2-obs2_climo;
error1_biased = mod1_biased-obs1_climo;
error2_biased = mod2_biased-obs2_climo;
[cost_mvn_uniform_before, cost_gmrf_uniform_before, cost_trad_uniform_before] ...
    = compute_cost( error1, error2, Q, alpha, sigma_mvn_obs, S_gmrf, p_q, figdir);
[cost_mvn_uniform_after, cost_gmrf_uniform_after, cost_trad_uniform_after] ...
    = compute_cost( error1_biased, error2_biased, Q, alpha, sigma_mvn_obs, S_gmrf, p_q, figdir);




close all

% What did we learn? Plot on same axis to compare types of errors. 

before_and_after_bias(cost_mvn_uniform_before, cost_gmrf_uniform_before, cost_trad_uniform_before,...
    cost_mvn_uniform_after, cost_gmrf_uniform_after, cost_trad_uniform_after)



summary_plots( ...
    cost_mvn_mean_perfect, cost_gmrf_mean_perfect, cost_trad_mean_perfect,...
    cost_mvn_mean_bias, cost_gmrf_mean_bias, cost_trad_mean_bias, ...
    cost_mvn_bias_xfield_cov, cost_gmrf_bias_xfield_cov, cost_trad_bias_xfield_cov,...
    cost_mvn_bias_infield_cov, cost_gmrf_bias_infield_cov, cost_trad_bias_infield_cov,... 
    cost_mvn_bias_inout_cov, cost_gmrf_bias_inout_cov, cost_trad_bias_inout_cov,...
    cost_mvn_bias_grad_one_field, cost_gmrf_bias_grad_one_field, cost_trad_bias_grad_one_field,...
    cost_mvn_bias_grad_both_field, cost_gmrf_bias_grad_both_field, cost_trad_bias_both_field,...
    cost_mvn_bias_grad_both_field_conf, cost_gmrf_bias_grad_both_field_conf, cost_trad_bias_grad_both_field_conf)


function [cost_mvn, cost_gmrf, cost_trad] = compute_cost(error1,error2, Q, alpha, sigma_mvn, S_gmrf, p_q, figdir)


B = eye( length( sigma_mvn ));
inv_sigma_mvn = sigma_mvn\B;
S_inv_gmrf  = S_gmrf\eye(2);
S_inv_trad  = diag(diag( S_gmrf)) \ eye(2);


% Reshape errors back to  points x time. 
T_error = reshape(error1, [],length(error1)); 
P_error = reshape(error2, [],length(error1)); 
errors = [T_error; P_error];

% modify Q using alpha 
Q_adj = alpha * eye(size(Q,1)) + (1- alpha)*Q; 
% The traditional and fields method use alpha=1
Q_trad = 1 * eye(size(Q,1)) + (1- 1)*Q;
Q_fields = Q_trad; 

% multiply times inverse spatially-avg covariance
% Don't need to do this for mvn because you already have the full sigma
% inverse. 
inv_sigma_gmrf = kron( S_inv_gmrf, Q_adj );
inv_sigma_trad = kron( S_inv_trad, Q_trad ); % trad uses alpha = 1 
inv_sigma_fields = kron( S_inv_gmrf, Q_fields ); % uses alpha = 1 but full S^-1 

% Compute each type of cost. 
% Now compute each experiment's errors and plot them using MVN. 
for i = 1:length(T_error)
  experiment = errors(:,i);
  cost_mvn(i)  = p_q(1) * experiment' * inv_sigma_mvn * experiment;
  cost_gmrf(i) = p_q(2) * experiment' * inv_sigma_gmrf * experiment;
  cost_trad(i) = p_q(3) * experiment' * inv_sigma_trad * experiment;
end

plots( error1,error2, cost_mvn,cost_gmrf, cost_trad, inv_sigma_mvn,inv_sigma_gmrf,inv_sigma_trad,inv_sigma_fields,figdir, sigma_mvn)

end


function plots(error1,error2,cost,cost_gmrf, cost_trad, inv_sigma_mvn,inv_sigma_gmrf,inv_sigma_trad,inv_sigma_fields,figdir, sigma_mvn )
mkdir(['figs/' , figdir ])


figure('Visible','off')
scatter(cost, cost_gmrf,'r'); hold on
scatter(cost, cost_trad,'b')
xlabel('cost MVN');ylabel('cost');
legend('GMRF','trad.')
outname = ['figs/' , figdir , '/cost_scatter.pdf'];
export_fig(outname)


fig=figure('Visible','off','Units','centimeters','Position',[1 1 20 17],'Color','w');
bottom = min( [inv_sigma_mvn(:)' inv_sigma_gmrf(:)' inv_sigma_trad(:)'] );
top = max( [inv_sigma_mvn(:)' inv_sigma_gmrf(:)' inv_sigma_trad(:)'] );
ax(1) = subplot(2,2,1); imagesc(inv_sigma_mvn); title( '\Sigma^{-1}' ); set(gca,'YTickLabel',[],'XTickLabel',[]); caxis([bottom top]);
ax(2) = subplot(2,2,2); imagesc( inv_sigma_gmrf ) ; title ( 'GMRF \Sigma^{-1}' ) ; set(gca,'YTickLabel',[],'XTickLabel',[]); caxis([bottom top])
ax(3) = subplot(2,2,3); imagesc( inv_sigma_fields ) ; title ( 'Fields \Sigma^{-1}' );set(gca,'YTickLabel',[],'XTickLabel',[]); caxis([bottom top])
ax(4) = subplot(2,2,4); imagesc( inv_sigma_trad ) ; title ( 'Trad. \Sigma^{-1}' );set(gca,'YTickLabel',[],'XTickLabel',[]); caxis([bottom top])
colormap(bluewhitered)
h=colorbar;
set(h, 'Position', [.86 .2 .03 .6])
for i=1:4
      pos=get(ax(i), 'Position');
      set(ax(i), 'Position', [0.85*pos(1) pos(2) pos(3) pos(4)]);
end
outname = ['figs/' , figdir , '/inv_sigma_pcolor.pdf'];
print(fig,'-dpdf', outname  )
%export_fig(outname)

fig=figure('Visible','off','Units','centimeters','Position',[1 1 12 15],'Color','w');
bottom = min( [inv_sigma_mvn(:)' inv_sigma_gmrf(:)' inv_sigma_trad(:)'] );
top = max( [inv_sigma_mvn(:)' inv_sigma_gmrf(:)' inv_sigma_trad(:)'] );
ax(1) = subplot(3,2,[1 2]); imagesc( sigma_mvn); title( '\Sigma' ); set(gca,'YTickLabel',[],'XTickLabel',[]); caxis([bottom top]);
ax(2) = subplot(3,2,3); imagesc(inv_sigma_mvn); title( '\Sigma^{-1}' ); set(gca,'YTickLabel',[],'XTickLabel',[]); caxis([bottom top]);
ax(3) = subplot(3,2,4); imagesc( inv_sigma_gmrf ) ; title ( 'GMRF \Sigma^{-1}' ) ; set(gca,'YTickLabel',[],'XTickLabel',[]); caxis([bottom top])
ax(4) = subplot(3,2,5); imagesc( inv_sigma_fields ) ; title ( 'Fields \Sigma^{-1}' );set(gca,'YTickLabel',[],'XTickLabel',[]); caxis([bottom top])
ax(5) = subplot(3,2,6); imagesc( inv_sigma_trad ) ; title ( 'Trad. \Sigma^{-1}' );set(gca,'YTickLabel',[],'XTickLabel',[]); caxis([bottom top])
colormap(bluewhitered)
h=colorbar('Location','south');
set(h, 'Position', [.1 .03 .8 .03])
for i=2:5
      pos=get(ax(i), 'Position');
      set(ax(i), 'Position', [pos(1) 1.05*pos(2) pos(3) pos(4)]);
end
%h.Label.String = 'Covariance';
%pos1=get(ax(1), 'Position');
set(ax(1),'Position',[0.35 0.70 0.35 0.25]); 
outname = ['figs/' , figdir , '/sigma_and_inv_sigma_pcolor.pdf'];
print(fig,'-dpdf', outname  )
export_fig(outname)



end

function plot_random_grid_point_corr( x1,x2, mod1,mod2)
% Validate that space and field correlations look right. 
% Choose a random grid cell in field 1 and plot correlations. 
i1 = randperm(length( x1 )); i2 = randperm(length( x2 )); 
ii=i1(1); jj=i2(1);
for i = 1:length(x1)
    for j = 1:length( x2)
        cf1(i,j) = corr2( squeeze( mod1(ii, jj,:)), squeeze(mod1(i,j,:)));
        cf2(i,j) = corr2( squeeze( mod1(ii, jj,:)), squeeze(mod2(i,j,:)));
    end
end

% Display this to see spatial and cross-field correlations to a random grid
% cell. Works as expected. 
figure('visible','off') 
colormap(bluewhitered)
subplot(1,2,1);
imagesc(cf1); 
h1=colorbar;
title(h1,'corr.')
caxis([0 1])
title('field1')
hold on
plot(jj,ii , 'r+', 'MarkerSize', 30)
subplot(1,2,2); 
imagesc(cf2); 
h2=colorbar;
title(h2,'corr.')
caxis([0  1])
title('field2')


% Other plots, commented for now

%figure
%ax1 = subplot(2,2,1); imagesc(mean( error1,3)); colormap( ax1, bluewhitered); colorbar; title('Ensemble Mean T error ')
%ax2 = subplot(2,2,2); imagesc(mean( error2,3)); colormap( ax2, bluewhitered); colorbar; title('Ensemble Mean P error ')
%ax3 = subplot(2,2,3); imagesc(error1(:,:,1));   colormap( ax3, bluewhitered); colorbar; title('T error member 1 ')
%ax4 = subplot(2,2,4); imagesc(error2(:,:,1));   colormap( ax4, bluewhitered); colorbar; title('P error member 1 ')
%outname = ['figs/' , figdir , '/syn_error.pdf'];
%export_fig(outname)

%figure
%max_cost = max( [max(cost) max(cost_gmrf) max(cost_trad) ] );
%bins = (0:10:max_cost);
%h1=histogram(cost_mvn,bins); hold on
%h2=histogram(cost_gmrf,bins); 
%h3=histogram(cost_trad,bins);
%ylabel('cost')
%legend('mvn','gmrf', 'trad')
%export_fig('cost_hist.png')


end

function summary_plots( ...
    cost_mvn_perfect, cost_gmrf_perfect, cost_trad_perfect,...
    cost_mvn_mean_bias, cost_gmrf_mean_bias, cost_trad_mean_bias, ...
    cost_mvn_bias_xfield_cov, cost_gmrf_bias_xfield_cov, cost_trad_bias_xfield_cov,...
    cost_mvn_bias_infield_cov, cost_gmrf_bias_infield_cov, cost_trad_bias_infield_cov,... 
    cost_mvn_bias_inout_cov, cost_gmrf_bias_inout_cov, cost_trad_bias_inout_cov,...
    cost_mvn_bias_grad_one_field, cost_gmrf_bias_grad_one_field, cost_trad_bias_grad_one_field,...
    cost_mvn_bias_grad_both_field, cost_gmrf_bias_grad_both_field, cost_trad_bias_both_field,...
    cost_mvn_bias_grad_both_field_conf, cost_gmrf_bias_grad_both_field_conf, cost_trad_bias_grad_both_field_conf)

sf1=figure('Units','centimeters','Position',[1 1 20 20]);
[M, T, G] = deal(cost_mvn_perfect,  cost_trad_perfect, cost_gmrf_perfect);
subplot(3,3,1);
scatter(M,T,'b');hold on
scatter(M,G,'r')
plot(linspace(0,max(M)), linspace(0, max(M)), 'k--')
xlabel('cost MVN')
ylabel('cost')
title('Perfect')
text(0.01, 0.6,  ['r =' num2str(round(corr2(M,T),2))], 'color','b', 'Units','normalized')
text(0.01, 0.5,  ['r =' num2str(round(corr2(M,G),2))], 'color','r', 'Units','normalized')
legend('trad.','GMRF','Location','northwest')
ylim([0 inf]); xlim([0 inf])
subplot(3,3,2);
[M, T, G] = deal(cost_mvn_mean_bias, cost_trad_mean_bias,cost_gmrf_mean_bias); 
scatter(M,T,'b');hold on
scatter(M,G,'r')
plot(linspace(0,max(M)), linspace(0, max(M)), 'k--')
xlabel('cost MVN')
ylabel('cost')
title('Mean bias')
text(0.01, 0.6,  ['r =' num2str(round(corr2(M,T),2))], 'color','b', 'Units','normalized')
text(0.01, 0.5,  ['r =' num2str(round(corr2(M,G),2))], 'color','r', 'Units','normalized')
legend('trad.','GMRF','Location','northwest')
ylim([0 inf]); xlim([0 inf])
subplot(3,3,3);
[M, T, G] = deal(cost_mvn_bias_xfield_cov, cost_trad_bias_xfield_cov,cost_gmrf_bias_xfield_cov); 
scatter(M,T ,'b');hold on
scatter(M,G,'r');
plot(linspace(0,max(M)), linspace(0, max(M)), 'k--')
xlabel('cost MVN')
ylabel('cost')
title('X-field correlation bias')
text(0.01, 0.6,  ['r =' num2str(round(corr2(M,T),2))], 'color','b', 'Units','normalized')
text(0.01, 0.5,  ['r =' num2str(round(corr2(M,G),2))], 'color','r', 'Units','normalized')
legend('trad.','GMRF','Location','northwest')
ylim([0 inf]); xlim([0 inf])
subplot(3,3,4);
[M, T, G] = deal(cost_mvn_bias_infield_cov, cost_trad_bias_infield_cov,cost_gmrf_bias_infield_cov ); 
scatter(M,T ,'b');hold on
scatter(M,G,'r');
plot(linspace(0,max(M)), linspace(0, max(M)), 'k--')
xlabel('cost MVN')
ylabel('cost')
title('In-field correlation bias')
text(0.01, 0.6,  ['r =' num2str(round(corr2(M,T),2))], 'color','b', 'Units','normalized')
text(0.01, 0.5,  ['r =' num2str(round(corr2(M,G),2))], 'color','r', 'Units','normalized')
legend('trad.','GMRF','Location','northwest')
ylim([0 inf]); xlim([0 inf])
subplot(3,3,5);
[M, T, G] = deal(cost_mvn_bias_inout_cov, cost_trad_bias_inout_cov,cost_gmrf_bias_inout_cov ); 
scatter(M,T ,'b');hold on
scatter(M,G,'r');
plot(linspace(0,max(M)), linspace(0, max(M)), 'k--')
xlabel('cost MVN')
ylabel('cost')
title('In- and x-field correlation bias')
text(0.01, 0.6,  ['r =' num2str(round(corr2(M,T),2))], 'color','b', 'Units','normalized')
text(0.01, 0.5,  ['r =' num2str(round(corr2(M,G),2))], 'color','r', 'Units','normalized')
legend('trad.','GMRF','Location','northwest')
ylim([0 inf]); xlim([0 inf])
subplot(3,3,6);
[M, T, G] = deal(cost_mvn_bias_grad_one_field, cost_trad_bias_grad_one_field,cost_gmrf_bias_grad_one_field ); 
scatter(M,T ,'b');hold on
scatter(M,G,'r');
plot(linspace(0,max(M)), linspace(0, max(M)), 'k--')
xlabel('cost MVN')
ylabel('cost')
title('Bias grad. 1 field')
text(0.01, 0.6,  ['r =' num2str(round(corr2(M,T),2))], 'color','b', 'Units','normalized')
text(0.01, 0.5,  ['r =' num2str(round(corr2(M,G),2))], 'color','r', 'Units','normalized')
legend('trad.','GMRF','Location','northwest')
ylim([0 inf]); xlim([0 inf])
subplot(3,3,7);
[M, T, G] = deal(cost_mvn_bias_grad_both_field, cost_trad_bias_both_field,cost_gmrf_bias_grad_both_field );
scatter(M,T ,'b');hold on
scatter(M,G,'r');
plot(linspace(0,max(M)), linspace(0, max(M)), 'k--')
xlabel('cost MVN')
ylabel('cost')
title('Consistent bias grad. 2 fields')
text(0.01, 0.6,  ['r =' num2str(round(corr2(M,T),2))], 'color','b', 'Units','normalized')
text(0.01, 0.5,  ['r =' num2str(round(corr2(M,G),2))], 'color','r', 'Units','normalized')
legend('trad.','GMRF','Location','northwest')
ylim([0 inf]); xlim([0 inf])
subplot(3,3,8);
[M, T, G] = deal( cost_mvn_bias_grad_both_field_conf, cost_trad_bias_grad_both_field_conf, cost_gmrf_bias_grad_both_field_conf);
scatter(M,T ,'b');hold on
scatter(M,G,'r');
plot(linspace(0,max(M)), linspace(0, max(M)), 'k--')
xlabel('cost MVN')
ylabel('cost')
title('Conflicting bias grad. 2 fields')
text(0.01, 0.6,  ['r =' num2str(round(corr2(M,T),2))], 'color','b', 'Units','normalized')
text(0.01, 0.5,  ['r =' num2str(round(corr2(M,G),2))], 'color','r', 'Units','normalized')
legend('trad.','GMRF','Location','northwest')
ylim([0 inf]); xlim([0 inf])

% re-plot 4 figs for paper: unbiased, mean bias, in-field bias, x-field bias, consistent
% bias two fields. These figs generated with:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%xfield_c = 0.5;                    % [0 to 0.9] cross-field correlation at the grid point
%infield_sc = 0.9;                   % [0 to 1] scale the spatial covariance within fields.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sf1=figure('Units','centimeters','Position',[1 1 15 15],'Color','w');
subplot(2,2,1);
[M, T, G] = deal(cost_mvn_perfect,  cost_trad_perfect, cost_gmrf_perfect); 
scatter(M,T,'b');hold on
scatter(M,G,'r')
plot(linspace(0,max(M)), linspace(0, max(M)), 'k--')
xlabel('E''(m) MVN')
ylabel('E''(m)')
title('Unbiased')
text(-0.2, 1, 'a', 'Units','Normalized','FontWeight','Bold','FontSize',16)
text(0.03, 0.6,  ['r = ' num2str(round(corr2(M,T),2))], 'color','b', 'Units','normalized')
text(0.03, 0.5,  ['r = ' num2str(round(corr2(M,G),2))], 'color','r', 'Units','normalized')
legend('Trad.','GMRF','Location','northwest')
ylim([0 inf]); xlim([0 inf])
subplot(2,2,2);
[M, T, G] = deal(cost_mvn_mean_bias, cost_trad_mean_bias,cost_gmrf_mean_bias); 
scatter(M,T,'b');hold on
scatter(M,G,'r')
plot(linspace(0,max(M)), linspace(0, max(M)), 'k--')
xlabel('E''(m) MVN')
ylabel('E''(m)')
title('Mean bias')
text(-0.2, 1, 'b', 'Units','Normalized','FontWeight','Bold','FontSize',16)
text(0.03, 0.6,  ['r = ' num2str(round(corr2(M,T),2))], 'color','b', 'Units','normalized')
text(0.03, 0.5,  ['r = ' num2str(round(corr2(M,G),2))], 'color','r', 'Units','normalized')
legend('Trad.','GMRF','Location','northwest')
ylim([0 inf]); xlim([0 inf])
subplot(2,2,3);
[M, T, G] = deal(cost_mvn_bias_xfield_cov, cost_trad_bias_xfield_cov,cost_gmrf_bias_xfield_cov);
scatter(M,T ,'b');hold on
scatter(M,G,'r');
plot(linspace(0,max(M)), linspace(0, max(M)), 'k--')
xlabel('E''(m) MVN')
ylabel('E''(m)')
title('Bias in x-field cov.')
text(-0.2, 1, 'c', 'Units','Normalized','FontWeight','Bold','FontSize',16)
text(0.03, 0.6,  ['r = ' num2str(round(corr2(M,T),2))], 'color','b', 'Units','normalized')
text(0.03, 0.5,  ['r = ' num2str(round(corr2(M,G),2))], 'color','r', 'Units','normalized')
legend('Trad.','GMRF','Location','northwest')
ylim([0 inf]); xlim([0 inf])
subplot(2,2,4);
[M, T, G] = deal(cost_mvn_bias_inout_cov, cost_trad_bias_inout_cov,cost_gmrf_bias_inout_cov ); 
scatter(M,T ,'b');hold on
scatter(M,G,'r');
plot(linspace(0,max(M)), linspace(0, max(M)), 'k--')
xlabel('E''(m) MVN')
ylabel('E''(m)')
title('Bias in x-field and in-field cov.')
text(-0.2, 1, 'd', 'Units','Normalized','FontWeight','Bold','FontSize',16)
text(0.03, 0.6,  ['r = ' num2str(round(corr2(M,T),2))], 'color','b', 'Units','normalized')
text(0.03, 0.5,  ['r = ' num2str(round(corr2(M,G),2))], 'color','r', 'Units','normalized')
legend('Trad.','GMRF','Location','northwest')
ylim([0 inf]); xlim([0 inf])
export_fig('figs/scatter_synthetics_four_cases.pdf')





figure('Units','centimeters','Position',[10 10 20 20])
ymax = 20;
subplot(3,3,1);
[T, G] = deal( cost_trad_perfect, cost_gmrf_perfect);
scatter(T,G)
text(0.6, 0.1,  ['r =' num2str(round(corr2(T,G),2))], 'Units','normalized')
xlabel('cost trad.')
ylabel('cost GMRF')
title('Perfect')
ylim([0 inf]); xlim([0 inf])
subplot(3,3,2);
[T, G] = deal( cost_trad_mean_bias ,cost_gmrf_mean_bias);
scatter(T,G)
text(0.6, 0.1,  ['r =' num2str(round(corr2(T,G),2))], 'Units','normalized')
xlabel('cost trad.')
ylabel('cost GMRF')
title('Mean bias')
ylim([0 inf]); xlim([0 inf])
subplot(3,3,3);
[T, G] = deal( cost_trad_bias_xfield_cov, cost_gmrf_bias_xfield_cov);
scatter(T,G)
text(0.6, 0.1,  ['r =' num2str(round(corr2(T,G),2))], 'Units','normalized')
xlabel('cost trad.')
ylabel('cost GMRF')
title('X-field correlation bias')
ylim([0 ymax])
ylim([0 inf]); xlim([0 inf])
subplot(3,3,4);
[T, G] = deal(cost_trad_bias_infield_cov,cost_gmrf_bias_infield_cov);
scatter(T,G)
text(0.6, 0.1,  ['r =' num2str(round(corr2(T,G),2))], 'Units','normalized')
xlabel('cost trad.')
ylabel('cost GMRF')
title('In-field correlation bias')
ylim([0 inf]); xlim([0 inf])
subplot(3,3,5);
[T, G] = deal(cost_trad_bias_inout_cov,cost_gmrf_bias_inout_cov);
scatter(T,G)
text(0.6, 0.1,  ['r =' num2str(round(corr2(T,G),2))], 'Units','normalized')
xlabel('cost trad.')
ylabel('cost GMRF')
title('IN- and x-field correlation bias')
ylim([0 inf]); xlim([0 inf])
subplot(3,3,6);
[T, G] = deal(cost_trad_bias_grad_one_field,cost_gmrf_bias_grad_one_field);
scatter(T,G)
text(0.6, 0.1,  ['r =' num2str(round(corr2(T,G),2))], 'Units','normalized')
xlabel('cost trad.')
ylabel('cost GMRF')
title('Bias grad. 1 field')
ylim([0 inf]); xlim([0 inf])
subplot(3,3,7);
[T, G] = deal(cost_trad_bias_both_field,cost_gmrf_bias_grad_both_field);
scatter(T,G)
text(0.6, 0.1,  ['r =' num2str(round(corr2(T,G),2))], 'Units','normalized')
xlabel('cost trad.')
ylabel('cost GMRF')
title('Consistent bias grad. 2 fields')
ylim([0 inf]); xlim([0 inf])
subplot(3,3,8);
[T, G] = deal(cost_trad_bias_grad_both_field_conf,cost_gmrf_bias_grad_both_field_conf);
scatter(T,G)
text(0.6, 0.1,  ['r =' num2str(round(corr2(T,G),2))], 'Units','normalized')
xlabel('cost trad.')
ylabel('cost GMRF')
title('Conflicting bias grad. 2 fields')
ylim([0 inf]); xlim([0 inf])
end

function before_and_after_bias( cost_mvn_uniform_before, cost_gmrf_uniform_before, cost_trad_uniform_before,...
    cost_mvn_uniform_after, cost_gmrf_uniform_after, cost_trad_uniform_after) 
% Is GMRF sensitive to mean bias? plot before and after. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%xfield_c = 0.5;                    % [0 to 0.9] cross-field correlation at the grid point
%infield_sc = 0.9;                   % [0 to 1] scale the spatial covariance within fields.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Units','centimeters','Position',[10 10 10.5 7],'Color','w');
[Mb,Tb,Gb] = deal(cost_mvn_uniform_before,  cost_trad_uniform_before, cost_gmrf_uniform_before);
[M, T, G] = deal(cost_mvn_uniform_after, cost_trad_uniform_after,cost_gmrf_uniform_after); 
scatter(Tb,T,'b');hold on
scatter(Gb,G,'r')
plot(linspace(0,max([Gb Tb])), linspace(0,max([Gb Tb])), 'k--')
xlabel('Unbiased mean')
ylabel('Biased mean')
title('Cost')
legend('Trad.','GMRF','Location','eastoutside')
ylim([0 inf]); xlim([0 inf])
export_fig('figs/before_and_after_bias.pdf')
end

