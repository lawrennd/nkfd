% DEMNKFD1 Demonstration of noise model with Gaussian distributions. First demo in ICML paper.
%
% COPYRIGHT : Neil D. Lawrence, 2000

% NKFD

rand('seed', 4e5);
randn('seed', 4e5);
%close all
%clear all
ndata = 200;
flipProb = 0.2;
dataMix = gmm(2, 2, 'full');
%dataMix.centres = [1 -1; -1 1; 1 1];
dataMix.centres = [1 -1; -2 1];
dataMix.covars = dataMix.covars;
C1 = [4 0; 2 4];
C2 = [2 0; 4 2];
dataMix.covars(:, :, 1) = C1*C1'
dataMix.covars(:, :, 2) = C2*C2'
dataMix.priors = [0.8 0.2];
[x, trueT] = gmmsamp(dataMix, ndata);
[x2, trueT2] = gmmsamp(dataMix, ndata);
trueT2 = trueT2>1;
trueT = trueT>1;

modelMix = gmm(2, 2, 'full');


noiseT = trueT;
noise = rand(size(trueT));
index = find(noise<flipProb);
noiseT(index) = ~trueT(index);

figure
plot(x(find(noiseT==0), 1), x(find(noiseT==0), 2), 'y.')
hold on
plot(x(find(noiseT==1), 1), x(find(noiseT==1), 2), 'wx')
hold on
xLim = get(gca, 'XLim');
yLim = get(gca, 'YLim');
zeroAxes(gca)
if exist('printPlot', 'var') & printPlot
  print -depsc ../tex/diagrams/demNkfd1_1.eps
  pos = get(gcf, 'paperposition');
  set(gcf, 'paperposition', [pos(1) pos(2) pos(3)/2 pos(4)/2]);
  print -dpng ../tex/diagrams/demNkfd1_1.png
  set(gcf, 'paperposition', [pos(1) pos(2) pos(3) pos(4)]);
end

figure
plot(x(find(trueT==0), 1), x(find(trueT==0), 2), 'y.')
hold on
plot(x(find(trueT==1), 1), x(find(trueT==1), 2), 'wx')
hold on
zeroAxes(gca)
if exist('printPlot', 'var') & printPlot
  print -depsc ../tex/diagrams/demNkfd1_2.eps
  pos = get(gcf, 'paperposition');
  set(gcf, 'paperposition', [pos(1) pos(2) pos(3)/2 pos(4)/2]);
  print -dpng ../tex/diagrams/demNkfd1_2.png
  set(gcf, 'paperposition', [pos(1) pos(2) pos(3) pos(4)]);
end
counter = 0;
for alpha1 = [0.0 0.3];
  alpha2 = alpha1;
  
  counter = counter + 1;
  modelMix.covars(:, :, 1) = cov(x(find(noiseT==0), :));
  modelMix.covars(:, :, 2) = cov(x(find(noiseT==1), :));
  modelMix.centres(1, :) = mean(x(find(noiseT==0), :));
  modelMix.centres(2, :) = mean(x(find(noiseT==1), :));
  modelMix.priors(1) = 1-alpha1;
  modelMix.priors(2) = 1-alpha2;
  options = foptions; 
  options(2) = 1e-6;
  options(3) = 1e-6;
  options(14) = 1000;
  options(1) = 1;
  modelMix = nlem(modelMix, x, noiseT, options);
  figure
  hold on
  if counter == 1
    title('(a)')
    fileName = 'demNkfd1_3';
  else
    title('(b)')
    fileName = 'demNkfd1_4';
  end
  for i = 1:modelMix.ncentres
    [v,d] = eig(modelMix.covars(:,:,i));
    for j = 1:2
      % Ensure that eigenvector has unit length
      v(:,j) = v(:,j)/norm(v(:,j));
      start=modelMix.centres(i,:)-sqrt(d(j,j))*(v(:,j)');
      endpt=modelMix.centres(i,:)+sqrt(d(j,j))*(v(:,j)');
      linex = [start(1) endpt(1)];
      liney = [start(2) endpt(2)];
      line(linex, liney, 'Color', 'r', 'LineWidth', 1)
    end
    % Plot ellipses of one standard deviation
    theta = 0:0.02:2*pi;
    newx = sqrt(d(1,1))*cos(theta);
    newy = sqrt(d(2,2))*sin(theta);
    % Rotate ellipse axes
    ellipse = (v*([newx; newy]))';
    % Adjust centre
    ellipse = ellipse + ones(length(theta), 1)*modelMix.centres(i,:);
    ellipse1 = plot(ellipse(:,1), ellipse(:,2), 'r-');
    set(ellipse1, 'LineWidth', 2)
  end
 
  for i = 1:dataMix.ncentres
    [v,d] = eig(dataMix.covars(:,:,i));
    for j = 1:2
      % Ensure that eigenvector has unit length
      v(:,j) = v(:,j)/norm(v(:,j));
      start=dataMix.centres(i,:)-sqrt(d(j,j))*(v(:,j)');
      endpt=dataMix.centres(i,:)+sqrt(d(j,j))*(v(:,j)');
      linex = [start(1) endpt(1)];
      liney = [start(2) endpt(2)];
      line(linex, liney, 'Color', 'r', 'LineWidth', 2, 'Linestyle', ':')
    end
    % Plot ellipses of one standard deviation
    theta = 0:0.02:2*pi;
    newx = sqrt(d(1,1))*cos(theta);
    newy = sqrt(d(2,2))*sin(theta);
    % Rotate ellipse axes
    ellipse = (v*([newx; newy]))';
    % Adjust centre
    ellipse = ellipse + ones(length(theta), 1)*dataMix.centres(i,:);
    ellipse2 = plot(ellipse(:,1), ellipse(:,2), 'r:');
    set(ellipse2, 'LineWidth', 2)
  end
  set(gca, 'XLim', xLim)
  set(gca, 'YLim', yLim)
  zeroAxes(gca)
  if exist('printPlot', 'var') & printPlot
    print('-depsc', ['../tex/diagrams/' fileName '.eps'])
    pos = get(gcf, 'paperposition');
    set(gcf, 'paperposition', [pos(1) pos(2) pos(3)/2 pos(4)/2]);
    print('-dpng', ['../tex/diagrams/' fileName '.png'])
    set(gcf, 'paperposition', [pos(1) pos(2) pos(3) pos(4)]);
  end
  hold off
  %/~
  % post = nlpost(modelMix, x, noiseT);
  
  %sum((post(:, 1) < .5) == trueT)
  
  % post = nlpost(modelMix, x2, rand(size(trueT2))>.5);  
  %sum((post(:, 1) < .5) == trueT2)
  %~/
end
