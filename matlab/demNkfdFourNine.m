% DEMNKFDFOURNINE Compare classification of USPS 4 vs 9 for a range of noise values.
%
% COPYRIGHT : Neil D. Lawrence, 2000

% NKFD

niters = 15;

classNum = 2;
[x, trueT, xTest, testT] = mapLoadData('fourNine');

% Make labels 0/1
trueT=trueT==1;
trueT = [trueT ~trueT];

testT=testT==1;
testT = [testT ~testT];
kern = kernCreate(x, 'rbf');
kern.inverseWidth = 1/128;
      
numData = size(x, 1);
numTestData = size(xTest, 1);
%/~
% x = [];
% numClass = zeros(classNum, 1);
% for i =1:classNum
%   X = [];
%   load([directory 'digit' num2str(dig(i))]);
%   %X = X(1:100, :);
%   x = [x; X];
%   numClass(i) = size(X, 1);
% end
  
% x = double(x);
% x = -(x-128)/128; 

% numData = sum(numClass);

% trueT = [];
% for i = 1:classNum
%   newT = zeros(numClass(i), classNum);
%   newT(:, i) = 1;
%   trueT = [trueT; newT];
% end
%~/
for alphaConst = [0.1 0]
  if alphaConst == 0
    numIters = 1;
  else 
    numIters = niters;
  end
  for seed = 1:10
    rand('seed',seed*1e5)
    randn('seed', seed*1e5)
    flipCount = 0;
    for flipProb = [0 .1 .2 .3 .4:0.02:.6 .7 .8 .9 1];
      flipCount = flipCount + 1;
      
      noiseT = trueT;
      for i = 1:numData
	if rand(1) < flipProb
	  temp = noiseT(i, :);
	  flipToZero = find(temp == 1);
	  index = find(temp == 0);
	  flipToOne = index(ceil(rand(1)*length(index)));
	  noiseT(i, flipToOne) = 1;
	  noiseT(i, flipToZero) = 0;
	end
      end
      
      numDataIn = 256;
      numTrainData = size(x, 1);
      alpha = ones(classNum, 1)*alphaConst;%.01;
      numKernelData = 200;
      [void, order] = sort(rand(1, numTrainData));
      kernelPoints = sort(order(1:numKernelData));
      xKern = x(kernelPoints, :);
      kx = kernCompute(kern, x, xKern);      
      hold on
      
      mix = gmm(numKernelData, classNum, 'full');
      
      
      mix.regularise = 0.0000000001;
      for i = 1:classNum
	switch mix.covar_type
	 case 'diag'
	  mix.covars(i, :) = diag(cov(kx(find(noiseT(:, i)==1), :))+mix.regularise* ...
				  eye(numKernelData));
	 case 'full'
	  mix.covars(:, :, i) = cov(kx(find(noiseT(:, i)==1), :))+mix.regularise* ...
	      eye(numKernelData);
	end
	mix.covars = mean(mix.covars, 3);
	mix.centres(i, :) = mean(kx(find(noiseT(:, i)==1), :));
	mix.priors(i) = 1;
      end
      options = foptions; 
      options(14) = 1;
      options(5) = 0;
      options(1) = 1;
      mix = nkfdEm(mix, kx, noiseT, options);
      
      for i = 1:numIters
	if i < 4
	  for i = 1:classNum
	    mix.priors(i) = 1-alpha(i);
	  end
	end
	mix = nkfdEm(mix, kx, noiseT, options);
	
	fprintf('Compare model against true class labels.\n')
	post = nkfdPost(mix, kx);
	predClass = (post>(1/mix.ncentres));
	fprintf('Train set %f per cent correctly labeled.\n',100* ...
		sum(noiseT(find(predClass==1)))/numData)
	
	switch mix.covar_type
	 case 'diag'
	  alphas = 1./mix.covars(1, :)'.*(mix.centres(1, :)' - mix.centres(2, :)');
	 case 'full'
	  alphas = inv(mix.covars(:, :, 1))*(mix.centres(1, :)' - mix.centres(2, :)');
	end
      end
      
      figure(1)
      clf
      width = 0;
      total = 0;
      for i = 1:classNum
	misClassed{i} = find(post(:, i)>1/classNum & noiseT(:, i)==1);
	numMis = length(misClassed{i});
	total = total + numMis;
	if numMis > width 
	  width = numMis;
	end
      end
      rowPlot = ceil(sqrt(total))+classNum+2;
      colPlot = ceil(sqrt(total));
      counter = 0;
    
      
      
      
      
      
      % Now compute the activations
      kxtest = kernCompute(kern, xTest, xKern); 
      testPost = nkfdPost(mix, kxtest);
      testPredClass = (testPost>(1/mix.ncentres));
      testErr(flipCount, seed) = 100*(sum(testT(find(testPredClass==1)))/numTestData);
      fprintf('Test set %f per cent correctly labeled.\n', testErr(flipCount, ...
						  seed))
    end
  end
  



  flipProb = [0 .1 0.2 .3 0.4:0.02:0.6 .7 0.8 .9 1];
  save(['uspsSigmoid9and4_alpha0_' ...
	num2str(alphaConst*10) '_noise.mat'],  'testErr', 'flipProb')

end

% Plot the rate of change of error as noise is increased
clf
load([directory 'uspsSigmoid9and4_noise.mat']);
linePlot = plot(flipProb, mean(testErr,2), '-')
set(linePlot, 'linewidth', 2);
load([directory 'uspsSigmoid9and4_no_noise.mat']);
hold on
linePlot = plot(flipProb, mean(testErr,2), ':')
set(linePlot, 'linewidth', 2);
ylabel('test % correct')
xlabel('label noise')
%plot(flipProb, binentropy(flipProb))
set(gca, 'Position', [0.13 0.41 0.775 0.415])
if exist('printPlot', 'var') & printPlot
  print -depsc ../tex/diagrams/demNkfdFourNine_1.eps
  pos = get(gcf, 'paperposition');
  set(gcf, 'paperposition', [pos(1) pos(2) pos(3)/2 pos(4)/2]);
  print -dpng ../tex/diagrams/demNkfdFourNine_1.png
  set(gcf, 'paperposition', [pos(1) pos(2) pos(3) pos(4)]);
end
figure
load([directory 'uspsSigmoid9and4_noise.mat']);
linePlot = plot(flipProb, mean(testErr,2), '-')
set(linePlot, 'linewidth', 2);
load([directory 'uspsSigmoid9and4_no_noise.mat']);
hold on
linePlot = plot(flipProb, mean(testErr,2), ':');
set(linePlot, 'linewidth', 2);
ylabel('test % correct')
xlabel('label noise')
%plot(flipProb, binentropy(flipProb))
set(gca, 'Position', [0.13 0.41 0.775 0.415])
set(gca, 'XLim', [0.4 0.6])
if exist('printPlot', 'var') & printPlot
  print -depsc ../tex/diagrams/demNkfdFourNine_2.eps
  pos = get(gcf, 'paperposition');
  set(gcf, 'paperposition', [pos(1) pos(2) pos(3)/2 pos(4)/2]);
  print -dpng ../tex/diagrams/demNkfdFourNine_2.png
  set(gcf, 'paperposition', [pos(1) pos(2) pos(3) pos(4)]);
end
