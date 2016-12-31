function gist = STgist(imgpath, param, gistpath)

param.G = createGabor(param.orientationsPerScale, param.imageSize+2*param.boundaryExtension);
% Precompute filter transfert functions (only need to do this once, unless
% image size is changes):
Nfeatures = size(param.G,3)*param.numberBlocks^2;

Imgs = imdir(imgpath);
Nscenes = length(Imgs);

% Loop: Compute gist features for all scenes
gist = zeros([Nscenes Nfeatures], 'single');
for n = 1:Nscenes
    g = [];   
    % if gist has already been computed, just read the file
    filegist = fullfile(gistpath, [Imgs(n).name(1:end-4) '.mat']);
    if exist(filegist, 'file')
        load(filegist, 'g');
    else
        % load image
        [img, ~] = STimread(imgpath,n);
        
        % convert to gray scale
        img = single(mean(img,3));

        % resize and crop image to make it square
        %img = imresizecrop(img, param.imageSize, 'bilinear');
        img = imresize(img, [param.imageSize, param.imageSize], 'bilinear'); %jhhays

        % scale intensities to be in the range [0 255]
        img = img-min(img(:));
        img = 255*img/max(img(:));
        
        % prefiltering: local contrast scaling
        output    = prefilt(img, param.fc_prefilt);

        % get gist:
        g = gistGabor(output, param);
        
        % save gist if a HOMEGIST file is provided

        save (filegist, 'g')

    end
    
    gist(n,:) = g;
    %drawnow
end


function output = prefilt(img, fc)
% ima = prefilt(img, fc);
% fc  = 4 (default)
% 
% Input images are double in the range [0, 255];
% You can also input a block of images [ncols nrows 3 Nimages]
%
% For color images, normalization is done by dividing by the local
% luminance variance.

if nargin == 1
    fc = 4; % 4 cycles/image
end

w = 5;
s1 = fc/sqrt(log(2));

% Pad images to reduce boundary artifacts
img = log(img+1);
img = padarray(img, [w w], 'symmetric');
[sn, sm, c, N] = size(img);
n = max([sn sm]);
n = n + mod(n,2);
img = padarray(img, [n-sn n-sm], 'symmetric','post');

% Filter
[fx, fy] = meshgrid(-n/2:n/2-1);
gf = fftshift(exp(-(fx.^2+fy.^2)/(s1^2)));
gf = repmat(gf, [1 1 c N]);

% Whitening
output = img - real(ifft2(fft2(img).*gf));
clear img

% Local contrast normalization
localstd = repmat(sqrt(abs(ifft2(fft2(mean(output,3).^2).*gf(:,:,1,:)))), [1 1 c 1]); 
output = output./(.2+localstd);

% Crop output to have same size than the input
output = output(w+1:sn-w, w+1:sm-w,:,:);



function g = gistGabor(img, param)
% 
% Input:
%   img = input image (it can be a block: [nrows, ncols, c, Nimages])
%   param.w = number of windows (w*w)
%   param.G = precomputed transfer functions
%
% Output:
%   g: are the global features = [Nfeatures Nimages], 
%                    Nfeatures = w*w*Nfilters*c

img = single(img);

w = param.numberBlocks;
G = param.G;
be = param.boundaryExtension;

if ndims(img)==2
    c = 1; 
    N = 1;
    [nrows ncols c] = size(img);
end
if ndims(img)==3
    [nrows ncols c] = size(img);
    N = c;
end
if ndims(img)==4
    [nrows ncols c N] = size(img);
    img = reshape(img, [nrows ncols c*N]);
    N = c*N;
end

[ny nx Nfilters] = size(G);
W = w*w;
g = zeros([W*Nfilters N]);

% pad image
img = padarray(img, [be be], 'symmetric');

img = single(fft2(img)); 
k=0;
for n = 1:Nfilters
    ig = abs(ifft2(img.*repmat(G(:,:,n), [1 1 N]))); 
    ig = ig(be+1:ny-be, be+1:nx-be, :);
    
    v = downN(ig, w);
    g(k+1:k+W,:) = reshape(v, [W N]);
    k = k + W;
    drawnow
end

if c == 3
    % If the input was a color image, then reshape 'g' so that one column
    % is one images output:
    g = reshape(g, [size(g,1)*3 size(g,2)/3]);
end


function y=downN(x, N)
% 
% averaging over non-overlapping square image blocks
%
% Input
%   x = [nrows ncols nchanels]
% Output
%   y = [N N nchanels]

nx = fix(linspace(0,size(x,1),N+1));
ny = fix(linspace(0,size(x,2),N+1));
y  = zeros(N, N, size(x,3));
for xx=1:N
  for yy=1:N
    v=mean(mean(x(nx(xx)+1:nx(xx+1), ny(yy)+1:ny(yy+1),:),1),2);
    y(xx,yy,:)=v(:);
  end
end

