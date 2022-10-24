
% Generate obs. 
% Input: x1, x2 = vectors for making meshgrid ( e.g. x1 = -2 :0.5: 2 )
%        ntime  = length of time series.  
%        cor    = the correlation coefficient in time between errors from different fields. 
function [f1_rshp, f2_rshp] = draw_from_mvn( x1, x2, mu, sigma_mvn, n )  
  
  random_vectors = mvnrnd( mu, sigma_mvn, n); % each row is a vector. 
  % uncomment if you want to see the correlations
  corrcoef(random_vectors);
  
  
  % reshape random vactors onto the grid.
  f1 = random_vectors(:,1:length(mu)/2);
  f2 = random_vectors(:,length(mu)/2 + 1:end);
  
  f1_rshp = reshape(f1', [length(x1),length(x1),n]);
  f2_rshp = reshape(f2', [length(x1),length(x1),n]);
   
 
end

