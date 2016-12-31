function [sift1] = ExtractSIFT(im1, pca_basis, sift_size)

% im = uint8(zeros(size(im1,1)+12,size(im1,2)+12,size(im1,3)));
% im(7:end-6,7:end-6,:) =im1;
% im(1:6,:,:) = repmat(im(7,:,:),[6,1,1]);
% im(end-5:end,:,:) = repmat(im(end-6,:,:),[6,1,1]);
% im(:,1:6,:) = repmat(im(:,6,:),[1,6,1]);
% im(:,end-5:end,:) = repmat(im(:,end-6,:),[1,6,1]);

im1 = imresize(im1, [size(im1,1)+12 size(im1,2)+12]); 
if ( nargin <= 2)
    sift_size = 4;
end

if  (ndims(im1) == 3 )
    im1 = rgb2gray(im1);
end
im1 = single(im1);


% Extract SIFT
[f1 sift1] = vl_dsift(im1, 'step', 1, 'size', sift_size, 'norm', 'FloatDescriptors', 'fast');

if ( ~isempty(pca_basis) ) % dimensionality reduction
    sift1 = pca_basis'*sift1;
end

% Formatting
x1 = unique(f1(1,:));
y1 = unique(f1(2,:));
width1 = numel(x1);
height1 = numel(y1);
sift1 = reshape(sift1', [height1, width1, size(sift1,1)]);
% mag1 = reshape(f1(3,:), [height1, width1]);

% lx = min(f1(1,:));
% rx = max(f1(1,:));
% ty = min(f1(2,:));
% dy = max(f1(2,:));
% sift = single(zeros(size(im1,1),size(im1,2),size(sift1,3)));
% sift(ty:dy,lx:rx,:) =sift1;
% sift(1:ty,:,:) = repmat(sift(ty,:,:),[ty,1,1]);
% sift(dy:end,:,:) = repmat(sift(dy,:,:),[size(im1,1)-dy+1,1,1]);
% sift(:,1:lx,:) = repmat(sift(:,lx,:),[1,lx,1]);
% sift(:,rx:end,:) = repmat(sift(:,rx,:),[1,size(im1,2)-rx+1,1]);
% bbox1 = [lx,rx,ty,dy]';
