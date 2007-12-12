function [post, a] = nkfdPost(mix, x, t)

% NKFDPOST Computes the class posterior probabilities of a KFD Model.
  
% NKFD

% Based on NETLAB code.
  

  
errstring = consist(mix, 'gmm', x);
if ~isempty(errstring)
  error(errstring);
end

ndata = size(x, 1);

[a, lna] = nkfdActiv(mix, x);

if nargin < 3
  priors = ones(ndata, mix.ncentres)/mix.ncentres;
else  
  
  priors = zeros(ndata, mix.ncentres);
  
  % Priors are associated with targets for noisy labels
  if size(t, 2) == 1
    index = find(t==1);
    if ~isempty(index);
      priors(index, :) = nrepmat([1-mix.priors(2) mix.priors(2)], 1, length(index));
    end
    index = find(t==0);
    if ~isempty(index)
      priors(index, :) = nrepmat([mix.priors(1) 1-mix.priors(1)], 1, length(index));
    end
  else
    for i = 1:mix.ncentres
      index = find(t(:, i)==1);
      priors(index, :) = (1-mix.priors(i))/(mix.ncentres-1);
      priors(index, i) = mix.priors(i);
    end
  end
end
priors(find(priors==0)) = eps;
lnpost = log(priors) + lna;
lnpost = lnpost - nrepmat(max(lnpost, [], 2), 2, mix.ncentres);
post = exp(lnpost);
%post = (ones(ndata, 1)*mix.priors).*a;
s = sum(post, 2);
% Set any zeros to one before dividing
s = s + (s==0);
post = post./(s*ones(1, mix.ncentres));







