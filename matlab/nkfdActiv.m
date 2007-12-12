function [a, lna] = nkfdActiv(mix, x)

% NKFDACTIV Computes the activations of a KFD Model.

% NKFD
  
% Based on Netlab code
  
% Check that inputs are consistent
errstring = consist(mix, 'gmm', x);
if ~isempty(errstring)
  error(errstring);
end

ndata = size(x, 1);

switch mix.covar_type

  case 'spherical'
    % Calculate squared norm matrix, of dimension (ndata, ncentres)
    n2 = dist2(x, mix.centres);

    % Calculate width factors
    wi2 = ones(ndata, 1) * (2 .* mix.covars);
    normal = (pi .* wi2) .^ (mix.nin/2);

    % Now compute the activations
    a = exp(-(n2./wi2))./ normal;

%/~
%  case 'diag'
%    a = zeros(ndata, mix.ncentres);
%    normal = (2*pi)^(mix.nin/2);
%    s = prod(sqrt(mix.covars), 2);
%    for i = 1:mix.ncentres
%      diffs = x - (ones(ndata, 1) * mix.centres(i, :));
%      a(:, i) = exp(-0.5*sum((diffs.*diffs)./(ones(ndata, 1) * ...
%	mix.covars(i,:)), 2)) ./ (normal*s(i));
%    end
%~/ 
 case 'diag'
    a = zeros(ndata, mix.ncentres);
    normal = mix.nin/2*log(2*pi);
    s = sum(.5*log(mix.covars), 2);
    for i = 1:mix.ncentres
      diffs = x - (ones(ndata, 1) * mix.centres(i, :));
      a(:, i) = -0.5*sum((diffs.*diffs)...
			 ./(ones(ndata, 1) * ...
			    mix.covars(i,:)), 2) - normal - s(i);
    end
    a = exp(a);
  case 'full'
    lna = zeros(ndata, mix.ncentres);
    lnnormal = mix.nin/2 * log(2*pi);
    [c, p] = chol(mix.covars(:,:,1));
    if p
      warning('Covariance is not positivie definite')
    end
    for i = 1:mix.ncentres
      diffs = x - (ones(ndata, 1) * mix.centres(i,:));
      % Use Cholesky decomposition of covariance matrix to speed computation
      temp = diffs/c;
      lna(:,i) = -0.5*sum(temp.*temp, 2) - lnnormal;
    end
    a = exp(lna);
end
  
