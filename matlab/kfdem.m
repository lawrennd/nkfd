function [mix, options, errlog] = kfdem(mix, x, t, options)

% KFDEM Kernel Fisher Discriminant EM algorithm.

% Based on NETLAB code.
  
% Check that inputs are consistent
errstring = consist(mix, 'gmm', x);
if ~isempty(errstring)
  error(errstring);
end

[ndata, xdim] = size(x);

% Sort out the options
if (options(14))
  niters = options(14);
else
  niters = 100;
end

display = options(1);
store = 0;
if (nargout > 2)
  store = 1;	% Store the error values to return them
  errlog = zeros(1, niters);
end
test = 0;
if options(3) > 0.0
  test = 1;	% Test log likelihood for termination
end

check_covars = 0;
if options(5) >= 1
  disp('check_covars is on');
  check_covars = 1;	% Ensure that covariances don't collapse
  MIN_COVAR = eps;	% Minimum singular value of covariance matrix
  init_covars = mix.covars;
end

% Main loop of algorithm
for n = 1:niters

  % Calculate posteriors based on old parameters
  [post, act] = kfdpost(mix, x, t);
  if size(t, 2) == 1
    index = find(t==1);
    post1 = [post(index, 2) post(index, 1)];
    priors(index, :) = nrepmat([1-mix.priors(2) mix.priors(2)], 1, length(index));
    index = find(t==0);
    post2 = post(index, :);
    priors(index, :) = nrepmat([mix.priors(1) 1-mix.priors(1)], 1, length(index));
  else
    for i = 1:mix.ncentres
      index = find(t(:, i) == 1);
      priors(index, :) = (1-mix.priors(i))/(mix.ncentres-1);
      priors(index, i) = mix.priors(i);
    end
  end
  
  %post = priors;%nrepmat(mix.priors, 1, ndata);
  %post2 = post(index, :);
  
  % Calculate error value if needed
  %  if (display | store | test)
  
  %    prob = sum(act.*priors, 2);
  %    % Error value is negative log likelihood of data
  %    index = find(prob<eps);
  %    prob(index) = eps;
  %    e = - sum(log(prob));
  %    if store
  %      errlog(n) = e;
  %    end
  %    if display > 0
  %      fprintf(1, 'Cycle %4d  Error %11.6f\n', n, e);
  %   end
  %    if test
  %      if (n > 1 & abs(e - eold) < options(3))
  %	options(8) = e;
  %	return;
  %      else
  %	eold = e;
  %      end
  %   end
  
  
  % Adjust the new estimates for the parameters
  
  new_pr = sum(post, 1);
  new_c = post' * x;
    
  % Now move new estimates to old parameter vectors
  if size(t, 2) == 1
    if mix.priors(1) ~= 1
      mix.priors(1) = mean(post2(:, 1));
    end
    if mix.priors(2) ~= 1
      mix.priors(2) = mean(post1(:, 1));
    end
  else
    for i = 1:mix.ncentres
      if mix.priors(i) ~= 1
	index = find(t(:, i) == 1);
	mix.priors(i) = mean(post(index, i));
      end
    end
  end
  disp(  mix.priors)
  %if mix.priors(1) < .5;
  %  mix.priors = 1-mix.priors;
  %  temp = mix.centres(1, :);
  %  mix.centres(1, :) = mix.centres(2, :);
  %  mix.centres(2, :) = temp;
  %  temp = mix.covars(:, :, 1);
  %  mix.covars(:, :, 1) = mix.covars(:, :, 2);
  %  mix.covars(:, :, 2) = temp;
  %  
  %end
  
  mix.centres = new_c ./ (new_pr' * ones(1, mix.nin));

  switch mix.covar_type
    case 'spherical'
      n2 = dist2(x, mix.centres);
      for j = 1:mix.ncentres
        v(j) = (post(:,j)'*n2(:,j));
      end
      mix.covars = ((v./new_pr))./mix.nin;
      if check_covars
	% Ensure that no covariance is too small
	for j = 1:mix.ncentres
	  if mix.covars(j) < MIN_COVAR
	    mix.covars(j) = init_covars(j);
	  end
	end
      end
    case 'diag'
     covar = zeros(size(mix.covars(1, :)));
     for j = 1:mix.ncentres
	diffs = x - (ones(ndata, 1) * mix.centres(j,:));
	covar = sum((diffs.*diffs).*(post(:,j)*ones(1, ...
	  mix.nin)), 1)./new_pr(j)+covar;
     end
     covar = covar/(ndata*mix.ncentres)+ mix.regularise*ones(1, mix.nin);
      if check_covars
	% Ensure that no covariance is too small
	for j = 1:mix.ncentres
	  if min(mix.covars(j,:)) < MIN_COVAR
	    mix.covars(j,:) = init_covars(j,:);
	  end
	end
      end
    case 'full'
     covar = zeros(mix.nin); 
     for j = 1:mix.ncentres       
       diffs = x - (ones(ndata, 1) * mix.centres(j,:));
       diffs = diffs.*(sqrt(post(:,j))*ones(1, mix.nin));
       %       covar = (diffs'*diffs)/new_pr(j) + covar;
       covar = (diffs'*diffs) + covar;
     end
     %     covar = covar/(mix.ncentres);
     covar = covar/(ndata*mix.ncentres)+ mix.regularise*eye(mix.nin);;
     for j = 1:mix.ncentres
       mix.covars(:, :, j) = covar;
     end
     if check_covars
       % Ensure that no covariance is too small
       for j = 1:mix.ncentres
	 if min(svd(mix.covars(:,:,j))) < MIN_COVAR
	   mix.covars(:,:,j) = init_covars(:,:,j);
	 end
       end
     end
  end

end

%options(8) = -sum(log(nlprob(mix, x, t)));
if (display >= 0)
  disp('Warning: Maximum number of iterations has been exceeded');
end

