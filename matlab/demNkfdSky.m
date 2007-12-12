% DEMKFDSKY Train a linear KFD with Label Noise on the Sky Data.

% NKFD

kernelType = 'linear';
kern = kernCreate(x, 'linear');
kernelParameter = 10;
close all
rand('seed', 2e5);
randn('seed', 2e5);

[xTrain, tTrain, xTest, tTest] = mapLoadData('sky');
ndata = size(xTrain, 1);
numInitIters = 1; % Number of iterations before mix.prior is adapted
numIters = 4; % Number of iterations adapting mix.prior

rand('seed', 4e5);
randn('seed', 4e5);


gamma1 = 0.0;
gamma2 = 0.5;



% How many points to use in kernel
numKernelData = 200;

% Select kernel points
[void, order] = sort(rand(1, ndata));
kernelPoints = sort(order(1:numKernelData));
xKern = xTrain(kernelPoints, :);

% Compute the kernel for training data
kx = kernCompute(kern, xTrain, xKern);
kk = kx(kernelPoints, :);

% initialise the kfd model
mix = gmm(numKernelData, 2, 'full');

% Set regularisation
mix.regularise = 1e-6;

% Set the covariance
covar1 = cov(kx(find(tTrain==0), :))+mix.regularise* ...
	 eye(numKernelData);
covar2 = cov(kx(find(tTrain==1), :))+mix.regularise* ...
	 eye(numKernelData);
covars = 0.5*(covar1 + covar2);
mix.covars(:, :, 1) = covars;
mix.covars(:, :, 2) = covars;

mix.centres(1, :) = mean(kx(find(tTrain==0), :));
mix.centres(2, :) = mean(kx(find(tTrain==1), :));

mix.priors(1) = 1 - gamma1;
mix.priors(2) = 1 - gamma2;

options = foptions; 
options(14) = 1;
options(5) = 0;
options(1) = 1;

for i = 1:(numIters + numInitIters)
  if i <= numInitIters
    mix.priors(1) = 1 - gamma1;
    mix.priors(2) = 1 - gamma2;
  end
  mix = nkfdEm(mix, kx, tTrain, options);
  
  %/~
  %fprintf('Test data against true class labels.\n')
  %post = nkfdPost(mix, kx);
  %fprintf('Train set %f per cent correctly labeled.\n',100* sum((post(:, 1) < .5) ==tTrain)/length(tTrain))
  %~/
end





covars = mix.covars(:, :, 1);

% Calculate the alphas for the Fisher discriminant
mns = mix.centres(1, :) - mix.centres(2, :);
cvs = mns'*mns;
[v1, d1] = eig(cvs, covars);
eigenVals = diag(d1);
[discEig, index] = max(eigenVals);
eigenVals(index) = [];
alphas = real(v1(:, index));
alphas = alphas/sqrt(alphas'*kk*alphas);
%/~
%S = covars; %cov(kx);
%outerAlphas = alphas*alphas';
%for i = 1:numPC
%  [v, d] = eig(S - kk*outerAlphas*S, kk);
%  eigVec = diag(d);
%  [maxEig, index] = max(eigVec);
%  eigVec(index) = [];
%  ndalphas{i} = real(v(:, index));
%  v(:, index) = [];
%  ndalphas{i} = ndalphas{i}/sqrt(ndalphas{i}'*kk*ndalphas{i});
%  outerAlphas = outerAlphas + ndalphas{i}*ndalphas{i}';
%end
%~/

% Load validation data
kxtest = kernCompute(kern, xTest, xKern);

testPost = nkfdPost(mix, kxtest);
testErr = sum((testPost(:, 1) < .5) == tTest)/length(tTest);
fprintf('Test set error %f per cent.\n',100*testErr)

skyclass('c:\datasets\construction\sky\training\63004.bmp', mix, 9, xKern, kernelType, kernelParameter);








