function nkfdClassVisualise(call)

% NKFDCLASSVISUALISE Sets up visualisation of decision boundary for NKFD model.
%
% COPYRIGHT : Neil D. Lawrence, 2000
  
% NKFD
  
global visualiseInfo
persistent called grayed

if isempty(called)
  grayed = 0;
end
called = 1;
%set(gcf, 'WindowButtonMotionFcn', 'classVisualise(''move'')')

switch call
 case 'click'
  oldK = visualiseInfo.k;
  visualiseInfo.k = 1;
  [x, y]  = localCheckPointPosition(visualiseInfo);  

  % If the pointer was not in the axes the values will be empty 
  if ~isempty(x)
    [dista, indexa] = localNearestPoint([x y], visualiseInfo, 1);
    [distb, indexb] = localNearestPoint([x y], visualiseInfo, 2);
    if dista < distb
      x = get(visualiseInfo.dataSet(1), 'XData');
      y = get(visualiseInfo.dataSet(1), 'YData');
      index = indexa;
      otherIndex = index;
    else
      x = get(visualiseInfo.dataSet(2), 'XData');
      y = get(visualiseInfo.dataSet(2), 'YData');
      index = indexb;
      otherIndex = index+visualiseInfo.numA(1);
    end
    x = x(index);
    y = y(index);
    %disp('Drawing')
    point = nkfdInvGetNormAxesPoint([x y], visualiseInfo.plotAxes);
    x = point(1);
    y = point(2);
    i = length(visualiseInfo.digitAxes);
    axisExist = otherIndex==visualiseInfo.digitIndex;
    if ~sum(axisExist) 
      visualiseInfo.digitAxes(i+1) = axes('position', ...
					  [x - visualiseInfo.axesWidth/2 ...
		    y - visualiseInfo.axesWidth/2 ...
		    visualiseInfo.axesWidth ...
		    visualiseInfo.axesWidth]);
      visualiseInfo.digitIndex(i+1) = otherIndex;
      image(reshape((visualiseInfo.x(otherIndex, :)+1)*32, 16, 16)', 'EraseMode', ...
	    'xor');
      set(visualiseInfo.digitAxes(i+1), 'visible', 'off')
      grayed = localGrayImage(visualiseInfo.digitAxes(i+1), otherIndex);
      colorMap gray
      %disp('drawing')
    else
      axisPos = find(axisExist);
      delete(visualiseInfo.digitAxes(axisPos));
      visualiseInfo.digitIndex(axisPos) = [];
      visualiseInfo.digitAxes(axisPos) = [];
      %disp('done already')
    end

    %    set(handle, 'visible', 'off')
    
  end
  visualiseInfo.k = oldK;
  
 case 'reinitialise'
  called = [];
  visualiseInfo = [];
  
 case 'editUpdated'
  newk = str2num(get(gcbo, 'string'));
  visualiseInfo = localUpdateK(visualiseInfo, newk); 
  
 case 'move'
  
  [x, y]  = localCheckPointPosition(visualiseInfo);  
  if ~isempty(x) 
    [dista, indexa] = localNearestPoint([x y], visualiseInfo, 1);
    [distb, indexb] = localNearestPoint([x y], visualiseInfo, 2);
    
    % numA is the number of points in class A
    indexb = indexb+visualiseInfo.numA(1);
    
    distAll = [dista; distb];
    [distAll, index] = sort(distAll);
    indexAll = [indexa; indexb];
    indexAll = indexAll(index(1:visualiseInfo.k));
    distAll = distAll(1:visualiseInfo.k);
    
    % Is pointer near an image?
    otherIndex = indexAll(1);
    if ~isempty(visualiseInfo.digitIndex)
      %disp('something is grayed')
      axisExist = otherIndex==visualiseInfo.digitIndex;
      axesPos = find(axisExist);
      % Is the nearby image not the old grayed image
      if sum(axisExist) 
	if otherIndex ~= grayed
	  if grayed ~= 0
	    % Ungrey the old image
	    otherAxesPos = find(grayed == visualiseInfo.digitIndex);
	    grayed = localUngrayImage(visualiseInfo.digitAxes(otherAxesPos));
	  end
	  % Grey the new one
	  grayed = localGrayImage(visualiseInfo.digitAxes(axesPos), otherIndex);
	end
      else
	if grayed ~= 0
	  otherAxesPos = find(grayed == visualiseInfo.digitIndex);
	  grayed = localUngrayImage(visualiseInfo.digitAxes(otherAxesPos));
	end
      end

    end
   
      % Invert distances
    dista = 1./dista;
    distb = 1./distb;
    distAll = 1./distAll;
    
    % Weighting factors for the images
    weightA = (dista./sum(dista))'*ones(1, 256);
    weightB = (distb./sum(distb))'*ones(1, 256);
    weightAll = (distAll./sum(distAll))'*ones(1, 256);
    
    % Draw the images of the nearest data points
    set(visualiseInfo.image(1), 'CData', reshape(sum(visualiseInfo.x(indexa, :).*weightA, 1), 16, 16)');
    set(visualiseInfo.image(2), 'CData', reshape(sum(visualiseInfo.x(indexb, :).*weightB, 1), 16, 16)');
    set(visualiseInfo.imageAll, 'CData', reshape(sum(visualiseInfo.x(indexAll, :).*weightAll, 1), 16, 16)');
  end
end

function visualiseInfo =localUpdateK(visualiseInfo, newk)

if visualiseInfo.k == newk
  return
elseif visualiseInfo.k < newk
  for i = visualiseInfo.k+1:newk
%    visualiseInfo.dataLine(1, i) = line([0 0], [0 0], 'color', [0 0 0]);
%    visualiseInfo.dataLine(2, i) = line([0 0], [0 0], 'color', [0 0 0]);
  end
elseif visualiseInfo.k > newk
  for i = visualiseInfo.k:-1:(newk+1)
%    set(visualiseInfo.dataLine(1, i), 'visible', 'off')
%    delete(visualiseInfo.dataLine(1, i))
%    set(visualiseInfo.dataLine(2, i), 'visible', 'off')
%    delete(visualiseInfo.dataLine(2, i))
  end
  
end
visualiseInfo.k = newk;



function [dist, index] = localNearestPoint(point, visualiseInfo, classLabel)

% Function to compute the nearest k points given current point, and a handle
% to data points (handle). lineHandle is a handle to the lines which are
% drawn to the nearest points.

x = get(visualiseInfo.dataSet(classLabel), 'xData')';
y = get(visualiseInfo.dataSet(classLabel), 'yData')';

points = [x y];
d = dist2(point, points);
[dist, index] = sort(d, 2);
index = index(1:visualiseInfo.k);
dist = dist(1:visualiseInfo.k);

for i = 1:visualiseInfo.k
%  set(visualiseInfo.dataLine(classLabel, i), 'XData', [point(1) x(index(i))]);
%  set(visualiseInfo.dataLine(classLabel, i), 'YData', [point(2) y(index(i))]);
end




function point = localGetNormCursorPoint(figHandle)

point = get(figHandle, 'currentPoint');
figPos = get(figHandle, 'Position');
% Normalise the point of the curstor
point(1) = point(1)/figPos(3);
point(2) = point(2)/figPos(4);

function [x, y] = localGetNormAxesPoint(point, axesHandle)

position = get(axesHandle, 'Position');
x = (point(1) - position(1))/position(3);
y = (point(2) - position(2))/position(4);
lim = get(axesHandle, 'XLim');
x = x*(lim(2) - lim(1));
x = x + lim(1);
lim = get(axesHandle, 'YLim');
y = y*(lim(2) - lim(1));
y = y + lim(1);


function grayed = localUngrayImage(axesHandle);

imageHandle = get(axesHandle, 'Children');
cdata = get(imageHandle, 'CData');
set(imageHandle, 'CData', (cdata-48)*4);
grayed = 0;

function grayed = localGrayImage(axesHandle, otherIndex);

imageHandle = get(axesHandle, 'Children');
cdata = get(imageHandle, 'CData');
set(imageHandle, 'CData', (cdata)/4+48);
grayed = otherIndex;

function [x, y] = localCheckPointPosition(visualiseInfo)

% Get the point of the cursor
point = localGetNormCursorPoint(gcf);

% get the position of the axes
position = get(visualiseInfo.plotAxes, 'Position');


% Check if the pointer is in the axes
if point(1) > position(1) ...
      & point(1) < position(1) + position(3) ...
      & point(2) > position(2) ...
      & point(2) < position(2) + position(4);
  
  % Rescale the point according to the axes
  [x y] = localGetNormAxesPoint(point, visualiseInfo.plotAxes);

  % Find the nearest point
else
  % Return nothing
  x = [];
  y = [];
end

