% Generate a covariance matrix and its inverse on a grid, for two
% correlated fields
% Input: x1, x2 = vectors for making meshgrid ( e.g. x1 = -2 :0.5: 2 )
          % doesnt actually use these, yet. assumes 3x3
%        cor    = the correlation coefficient in time between errors from different fields, at the same gridpoint. 
%        multiplier = a parameter for scaling the correlation within
%        fields. Set this too high, and the matrix becomes singular. 
function [sigma_mvn, inv_sigma_mvn] = make_sigma( x1, x2,xfield_c, infield_sc ) 
  % Generate the covariance matrix.
  % This is actually correlation, assuming standard normal vars. 
  % Rules are: border = 0.8, one cell between = 0.4, corners touch = 0.3
  
    in_field_bottom = ...
               [1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
                0.8 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
                0.4 0.8 1.0 0.0 0.0 0.0 0.0 0.0 0.0;
                0.8 0.3 0.0 1.0 0.0 0.0 0.0 0.0 0.0;
                0.3 0.8 0.3 0.8 1.0 0.0 0.0 0.0 0.0;
                0.0 0.3 0.8 0.3 0.8 1.0 0.0 0.0 0.0;
                0.4 0.0 0.0 0.8 0.3 0.0 1.0 0.0 0.0; 
                0.0 0.4 0.0 0.3 0.8 0.3 0.8 1.0 0.0;
                0.0 0.0 0.4 0.0 0.3 0.8 0.4 0.8 1.0];      
            
 
  % This original version is not positive semi-definite. Adjust
  % correlations downward
  in_field_bottom(in_field_bottom == 0.8) = 0.55;
  in_field_bottom(in_field_bottom == 0.4) = 0.15;
  in_field_bottom(in_field_bottom == 0.3) = 0.35;          
  % symmetry. 
  in_field = in_field_bottom + tril(in_field_bottom,-1)';
  % scale the in_field covariance 
  % Multiply by the infield multiplier. This scales the in-field covariance, so you only want it to apply there.  
  old_in_field_diags = diag(diag(eye(size( in_field))));
  in_field   = in_field * infield_sc ;
  % It's now invertible for negatives, but apparently not positive
  % semidefinite
  % keep the original diagonals. 
  new_in_field = in_field - diag(diag(in_field)) + old_in_field_diags;
   % add second field, correlated to the first field by the parameter "cor"
  cross_field = new_in_field * xfield_c;

  
  sigma_mvn = [new_in_field cross_field;  cross_field new_in_field];

  
  figure('visible','off','Color','w')
  imagesc( sigma_mvn)
  set(gca,'YTickLabel',[],'XTickLabel',[])
  colormap(bluewhitered);
  title('\Sigma')
  a = colorbar;
  a.Label.String = 'Covariance';
  export_fig('figs/covmat.png')
    
  % Lets make sure this covariance matrix can be inverted. 
  B = eye( length( sigma_mvn ));
  %inv_sigma_mvn = sigma_mvn\B;  
  %imagesc( inv_sigma_mvn); colorbar % uncomment to glimpse the inv cov matrix. 
  inv_sigma_mvn = inv( sigma_mvn);
 
  
  % Curious: how many positive eigenvalues does it have? 
  [V,D] = eig(sigma_mvn); % looks like all 18 are positive. 
