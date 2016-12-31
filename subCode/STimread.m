function [img ImgName]= STimread(imgpath,index)
Imgs = imdir(imgpath);
[~, ImgName] = fileparts(Imgs(index).name);
if exist(fullfile(imgpath, [ImgName '.png']),'file')           
	img = imread(fullfile(imgpath, [ImgName '.png']));
elseif exist(fullfile(imgpath, [ImgName '.jpg']),'file')
    img = imread(fullfile(imgpath, [ImgName '.jpg']));
elseif exist(fullfile(imgpath, [ImgName '.bmp']),'file')
    img = imread(fullfile(imgpath, [ImgName '.bmp']));
end 



