function [post, a] = nlpost(mix, x, t)

% NLPOST Computes the posteriors for a noiseylabel model

% NKFD
  
% Based on netlab code.
  
% Check that inputs are consistent
errstring = consist(mix, 'gmm', x);
if ~isempty(errstring)
  error(errstring);
end

ndata = size(x, 1);

a = gmmactiv(mix, x);

priors = zeros(ndata, mix.ncentres);

% Priors are associated with targets for noisey labels
if nargin > 2
  index = find(t==1);
  priors(index, :) = nrepmat([1-mix.priors(2) mix.priors(2)], 1, length(index));
  
  index = find(t==0);
  priors(index, :) = nrepmat([mix.priors(1) 1-mix.priors(1)], 1, length(index));
else
  priors = nrepmat([0.5 0.5], 1, size(x, 1));
end
post = priors.*a;
%post = (ones(ndata, 1)*mix.priors).*a;
s = sum(post, 2);
% Set any zeros to one before dividing
s = s + (s==0);
post = post./(s*ones(1, mix.ncentres));






