function saliency= saliencycompact(Sal,ratio)
[height,width] = size(Sal);
% hfactor = 0.4*height;
% wfactor = 0.4*width;
[X,Y] = meshgrid(1:width,1:height); 
T = exp(Sal);
X = X.*T;
Y = Y.*T;

salCenter = [sum(X(:))/sum(T(:)) sum(Y(:))/sum(T(:))];
[X,Y] = meshgrid(1:width,1:height);
X = (X-salCenter(1)).^2.*T;
Y = (Y-salCenter(2)).^2.*T;
avwidth = sum(X(:))/sum(T(:));
avheight = sum(Y(:))/sum(T(:));
avwidth = ratio*avwidth;
avheight = ratio*avheight;
% avwidth = exp(1-avwidth/wfactor)*avwidth;
% avheight = exp(1-avheight/hfactor)*avheight;
[X,Y] = meshgrid(1:width,1:height);
saliency = (exp(-((Y-salCenter(:,2)).^2/avheight+ (X-salCenter(:,1)).^2/avwidth)));

end