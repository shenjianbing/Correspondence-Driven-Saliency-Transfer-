%%%%%%%%%%%%%%ÒÀÀµ±³¾°¹À¼Æ+geo%%%%%%%%%%%%%%%
function [BSsaliency mask leftx rightx lefty righty]= centerBias(IniSal_region,IniSal,superpixels) 
[height,width] = size(IniSal);
nLabel = double(max(superpixels.Label(:)));
Label = superpixels.Label;
mask = (IniSal>mean(IniSal_region(:)));
mask = int32(mask>0);
fd = mask.*Label;
fd = unique(fd(:));
fd(fd==0) = [];
fd = sort(fd);
bd = int32(IniSal<=mean(IniSal_region(:))).*Label;
bd = unique(bd(:));
bd(bd==0) = [];
bd_f = IniSal_region(bd);

colDistM = squareform(pdist(superpixels.Lab));
colDistMx = colDistM;
[conSPix Conedge]= find_connect_superpixel( Label, nLabel, height ,width );  
ConSPix=sparse([Conedge(:,1);Conedge(:,2)],[Conedge(:,2);Conedge(:,1)], ...
         [ones(size(Conedge(:,1)));ones(size(Conedge(:,1)))],nLabel,nLabel);
ConSPix = full(ConSPix);
ConSPix = ConSPix +eye(size(ConSPix));
              
b_geoDist = GeodesicSaliency(ConSPix, double(bd), colDistM, 0,bd_f,[]);
R = b_geoDist';
R = (R-min(R))/(max(R)-min(R));
fp = reshape(superpixel2pixel(double(superpixels.Label),R),height,width);
x =12;







