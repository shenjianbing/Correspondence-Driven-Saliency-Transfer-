function Im=im2feature(im,ratio)
if exist('ratio')~=1
    ratio=1;
end
if isfloat(im)~=1
    im=im2double(im);
end
dxfilter=[1 -8 0 8 -1]/12;
dyfilter=dxfilter';
if size(im,3)==1
    Im=im*ratio;
%     Im(:,:,2)=imlocalenhance(imfilter(im,dxfilter,'same','replicate'));
%     Im(:,:,3)=imlocalenhance(imfilter(im,dyfilter,'same','replicate'));
    im=imfilter(im,fspecial('gaussian',5,1),'same','replicate');
    foo1=imfilter(im,dxfilter,'same','replicate');
    foo2=imfilter(im,dyfilter,'same','replicate');
    Im(:,:,2)=foo1;
    Im(:,:,3)=foo2;
%     H=sqrt(imfilter(foo1.^2+foo2.^2,fspecial('gaussian',15,4),'same','replicate'));
%     Im(:,:,2)=foo1./(H+0.01);
%     Im(:,:,3)=foo2./(H+0.01);
else
    Im=rgb2gray(im);
    Im(:,:,2)=imfilter(Im,dxfilter,'same','replicate');
    Im(:,:,3)=imfilter(Im(:,:,1),dyfilter,'same','replicate');
    Im(:,:,4)=(im(:,:,2)-im(:,:,1))*0.25;
    Im(:,:,5)=(im(:,:,2)-im(:,:,3))*0.25;
end

