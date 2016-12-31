%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sMap = morphSmooth(I,width)
% opening by reconstruction followed by closing by reconstruction
% see the following material for detailed explanations 
% http://www.mathworks.com/products/demos/image/watershed/ipexwatershed.html
I = uint8(I*255);
se = strel('square',width);
Ie = imerode(I, se);
Iobr = imreconstruct(Ie, I);
Iobrd = imdilate(Iobr, se);
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
sMap = mat2gray(Iobrcbr);