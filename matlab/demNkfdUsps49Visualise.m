% DEMNKFDUSPS49VISUALISE Visualise the decision boundary between USPS fours and nines.

global visualiseInfo
kernelType = 'rbf';
kernelParameter = 128;
close all
rand('seed', 2e5);
randn('seed', 2e5);

%[x, trueT] = mapLoadData('cedar69');
[x, trueT, xtest, ttest] = mapLoadData('fourNine');
trueT = trueT==1;
ttest = ttest==1;
numDataIn = size(x, 2);
numTrainData = size(x, 1);

flipProb = 0;

% Add noise to labels
noiseT = trueT;
noise = rand(size(trueT));
index = find(noise<flipProb);
noiseT(index) = ~trueT(index);

numPC = 2;
display = 0;
numInitIters = 0; % Number of iterations before mix.prior is adapted
numIters = 1; % Number of iterations adapting mix.prior
gamma = 0.1;

% Set prior flip probabilities
gamma1 = gamma;
gamma2 = gamma;


% How many points to use in kernel
numKernelData = 200;

% Select kernel points
[void, order] = sort(rand(1, numTrainData));
kernelPoints = sort(order(1:numKernelData));
xKern = x(kernelPoints, :);

% Compute the kernel for training data
kern = kernCreate(x, kernelType);
kern.inverseWidth = 1/128;
kx = kernCompute(kern, x, xKern);
%kernel(x, xKern, kernelType, kernelParameter);
kk = kx(kernelPoints, :);

% initialise the kfd model
mix = gmm(numKernelData, 2, 'full');
% Set regularisation
mix.regularise = 0.0000000001;

% Set the covariance
switch mix.covar_type
 case 'diag'
  covar1 = diag(cov(kx(find(noiseT==0), :))+mix.regularise* ...
			  eye(numKernelData));
  covar2 = diag(cov(kx(find(noiseT==1), :))+mix.regularise* ...
			  eye(numKernelData));
  covars = 0.5*(covar1 + covar2);
  mix.covars(1, :) = covars;
  mix.covars(2, :) = covars;
  
 case 'full'
  covar1 = cov(kx(find(noiseT==0), :))+mix.regularise* ...
      eye(numKernelData);
  covar2 = cov(kx(find(noiseT==1), :))+mix.regularise* ...
      eye(numKernelData);
  covars = 0.5*(covar1 + covar2);
  mix.covars(:, :, 1) = covars;
  mix.covars(:, :, 2) = covars;
end

mix.centres(1, :) = mean(kx(find(noiseT==0), :));
mix.centres(2, :) = mean(kx(find(noiseT==1), :));

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
  mix = nkfdEm(mix, kx, noiseT, options);
  
  fprintf('Test data against true class labels.\n')
  post = nkfdPost(mix, kx);
  fprintf('Train set %f per cent correctly labeled.\n',100* sum((post(:, 1) < .5) ==noiseT)/length(noiseT))
  
end

numTestData = size(xtest,1);

% Compute test kernel
kxtest = kernCompute(kern, xtest, xKern);

testPost = nkfdPost(mix, kxtest);
testErr = sum((testPost(:, 1) < .5) == ttest)/length(ttest);
fprintf('Test set error %f per cent.\n',100*testErr)

%/~
%switch mix.covar_type
% case 'diag'
%  alphas = 1./mix.covars(1, :)'.*(mix.centres(1, :)' - mix.centres(2, :)');
% case 'full'
%  alphas = inv(mix.covars(:, :, 1))*(mix.centres(1, :)' - mix.centres(2, :)');
%end
%~/
if display
  figure(1)
  clf
  misclassedOne = find(post(:, 1)>.5 & noiseT==1);
  misclassedTwo = find(post(:, 2)>.5 & noiseT==0);
  total = length(misclassedOne) + length(misclassedTwo);
  rowPlot = ceil(sqrt(total))+2;
  colPlot = ceil(sqrt(total));
  width = max([length(misclassedOne), length(misclassedTwo)]);
  counter = 0;
  
  for i = misclassedOne'
    counter = counter + 1;
    subplot(rowPlot, colPlot, counter)
    imagesc(reshape(x(i, :), 16, 16))
    axis image
    axis off
    colormap('gray')
    drawnow
  end
  counter = (1 +ceil(counter/colPlot))*colPlot
  for i = misclassedTwo'
    counter = counter + 1;
    subplot(rowPlot, colPlot, counter)
    imagesc(reshape(x(i, :), 16, 16))
    axis image
    axis off
    colormap('gray')
    drawnow
  end  
  
  figure(2)
  clf
  misclassedOne = find(testPost(:, 1)>.5 & ttest==1);
  misclassedTwo = find(testPost(:, 2)>.5 & ttest==0);
  total = length(misclassedOne) + length(misclassedTwo);
  rowPlot = ceil(sqrt(total))+2;
  colPlot = ceil(sqrt(total));
  width = max([length(misclassedOne), length(misclassedTwo)]);
  counter = 0;
  
  for i = misclassedOne'
    counter = counter + 1;
    subplot(rowPlot, colPlot, counter)
    imagesc(reshape(xtest(i, :), 16, 16))
    axis off
    axis image
    colormap('gray')
    drawnow
  end
  counter = (1 +ceil(counter/colPlot))*colPlot
  for i = misclassedTwo'
    counter = counter + 1;
    subplot(rowPlot, colPlot, counter)
    imagesc(reshape(xtest(i, :), 16, 16))
    axis off
    axis image
    colormap('gray')
    drawnow
  end
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

S = covars; %cov(kx);
outerAlphas = alphas*alphas';
for i = 1:numPC
  [v, d] = eig(S - kk*outerAlphas*S, kk);
  eigVec = diag(d);
  [maxEig, index] = max(eigVec);
  eigVec(index) = [];
  ndalphas{i} = real(v(:, index));
  v(:, index) = [];
  ndalphas{i} = ndalphas{i}/sqrt(ndalphas{i}'*kk*ndalphas{i});
  outerAlphas = outerAlphas + ndalphas{i}*ndalphas{i}';
end



% Calculate the alphas for first two kernel Principal Components
S = cov(kx);
[v, d] = eig(S, kk);
eigVec = diag(d);
for i = 1:numPC
  [maxEig, index] = max(eigVec);
  eigVec(index) = [];
  pcalphas{i} = real(v(:, index));
  v(:, index) = [];
  pcalphas{i} = pcalphas{i}/sqrt(pcalphas{i}'*kk*pcalphas{i});
end

[v, d] = eig(covars, kk);
eigVec = diag(d);
for i = 1:numPC
  [maxEig, index] = max(eigVec);
  eigVec(index) = [];
  fdalphas{i} = real(v(:, index));
  v(:, index) = [];
  fdalphas{i} = fdalphas{i}/sqrt(fdalphas{i}'*kk*fdalphas{i});
end


% Calculate the kernel for the first class
kx0 = kx(find(trueT==0), :);

% Calculate the kernel for the second class
kx1 = kx(find(trueT==1), :);


% Number of data in first class
visualiseInfo.numA(1) = sum((trueT==0));
visualiseInfo.axesWidth = 0.015;

% Create the image sub plots
figure(2)
imageAxesa = subplot(1, 3, 1);
visualiseInfo.image(1) = imagesc(reshape(x(1, :), 16, 16));
imageAxesb = subplot(1, 3, 3);
visualiseInfo.image(2) = imagesc(reshape(x(1+visualiseInfo.numA, :), 16, 16));
imageAxesc = subplot(1, 3, 2);
visualiseInfo.imageAll = imagesc(reshape(x(1+visualiseInfo.numA, :), 16, 16));
colormap gray
figure(1)
% Create the plot for the data
visualiseInfo.plotAxes = subplot(1, 1, 1);

% Determine projections of the data points
% X direction from Fisher Discriminant
x0 = kx0*alphas;
x1 = kx1*alphas;
% Y direction from uncorrelated Principal component
y0 = kx0*ndalphas{1};
y1 = kx1*ndalphas{1};

% Plot the data
visualiseInfo.dataSet(1) = plot(x0, y0, 'r.');
set(visualiseInfo.dataSet(1), 'MarkerSize', 10)
hold on
visualiseInfo.dataSet(2) = plot(x1, y1, 'bx');
set(visualiseInfo.dataSet(2), 'MarkerSize', 3)



% Draw the initial lines
visualiseInfo.k = 1; % Number of points to average over
%for i = 1:visualiseInfo.k
%  visualiseInfo.dataLine(1, i)= line([0 0], [1 1], 'EraseMode', 'xor');
%  visualiseInfo.dataLine(2, i)= line([0 0], [1 1], 'EraseMode', 'xor');
%end

% Set up the X limits and Y limits of the main plot
xLim = [min([x0; x1]) max([x0; x1])];
xSpan = xLim(2) - xLim(1);
xLim(1) = xLim(1) - 0.05*xSpan;
xLim(2) = xLim(2) + 0.05*xSpan;
xSpan = xLim(2) - xLim(1);

yLim = [min([y0; y1]) max([y0; y1])];
ySpan = yLim(2) - yLim(1);
yLim(1) = yLim(1) - 0.05*ySpan;
yLim(2) = yLim(2) + 0.05*ySpan;
ySpan = yLim(2) - yLim(1);

set(visualiseInfo.plotAxes, 'XLim', xLim)
set(visualiseInfo.plotAxes, 'YLim', yLim)
visualiseInfo.digitAxes = [];
visualiseInfo.digitIndex = [];

% Turn off the image axes
set(imageAxesa, 'visible', 'off')
set(imageAxesb, 'visible', 'off')
set(imageAxesc, 'visible', 'off')

% Pass the data to visualiseInfo
visualiseInfo.x = x;

hold off
figure(2)
editBox = uicontrol('Callback','nkfdClassVisualise(''editUpdated'')', ...
	'Position',[309 29.25 39 20.25], ...
	'String',num2str(visualiseInfo.k), ...
	'Style','edit', ...
	'Tag','kValue');
% Set the callback function
figure(1)
set(gcf, 'WindowButtonMotionFcn', 'nkfdClassVisualise(''move'')')
set(gcf, 'WindowButtonUpFcn', 'nkfdClassVisualise(''click'')')

