function [BSsaliency, mask]= centerBias_geodesic(IniSal_region,IniSal,superpixels)
ratio = 1;
[height,width] = size(IniSal);
nLabel = double(max(superpixels.Label(:)));
Label = superpixels.Label;
L{1} = uint32(Label);
bdIds=[Label(1,:)';Label(end,:)';Label(:,1);Label(:,end)];
bdIds = unique(bdIds);

if mean(IniSal_region(Label(1,:)'))>0.2
    IniSal_region(Label(1,:)') = IniSal_region(Label(1,:)')*0.5;
end
if mean(IniSal_region(Label(end,:)'))>0.2
    IniSal_region(Label(end,:)') = IniSal_region(Label(end,:)')*0.5;
end
if mean(IniSal_region(Label(:,1)))>0.2
    IniSal_region(Label(:,1)) = IniSal_region(Label(:,1))*0.5;
end
if mean(IniSal_region(Label(1,:)'))>0.2
    IniSal_region(Label(:,1)) = IniSal_region(Label(:,1))*0.5;
end

colDistM = squareform(pdist(superpixels.Lab));

[~, Conedge]= find_connect_superpixel2( Label, nLabel, height ,width );  
ConSPix=sparse([Conedge(:,1);Conedge(:,2)],[Conedge(:,2);Conedge(:,1)], ...
         [ones(size(Conedge(:,1)));ones(size(Conedge(:,1)))],nLabel,nLabel);
%ConSPix = full(ConSPix);
ConSPix = ConSPix -diag(diag(ConSPix));
clipVal = EstimateDynamicParas(ConSPix, colDistM);
ConSPix = ConSPix+eye(size(ConSPix));
              
b_geoDist = GeodesicSaliency(ConSPix, double(bdIds), colDistM, clipVal,IniSal_region(bdIds));
R = b_geoDist';%.*(0.5+IniSal_region);
R = (R-min(R))/(max(R)-min(R));
BSsaliency = reshape(superpixel2pixel(double(Label),R),height,width);
BSsaliency = (0.0+BSsaliency).*saliencycompact(BSsaliency,ratio);
S{1} = repmat(BSsaliency(:),[1 3]);
[ R, ~, ~ ] = getSuperpixelStats(S(1:1),L, nLabel);
R = double(R(:,1));
R = (R-min(R))/(max(R)-min(R));
BSsaliency = reshape(superpixel2pixel(double(Label),R),height,width);

mask = int32(BSsaliency>mean(BSsaliency(:)));
fd = mask.*Label;
fd = unique(fd(:));
fd(fd==0) = [];
fd = sort(fd);
fcenters = superpixels.centres( fd,:);
mask = zeros(height,width);
   %  K-means Algorithm 
    paramPropagate.nclus = 3;    %the number of boundary clusters
    paramPropagate.maxIter=200;
    centers_fg = form_codebook(fcenters',...
                 paramPropagate.nclus,paramPropagate.maxIter);
    [ featLabel_fg ] = labelCluster( centers_fg, ...
                 fcenters', length(fd), paramPropagate.nclus );    
    clus_num=numel(unique(featLabel_fg));
    maxlabel=max(featLabel_fg);
    if (clus_num~= paramPropagate.nclus)
        if clus_num == maxlabel 
            paramPropagate.nclus = clus_num;
        else
            paramPropagate.nclus = clus_num;
            centers_fg = form_codebook(fcenters',...
                         paramPropagate.nclus,paramPropagate.maxIter);
            [ featLabel_fg ] = labelCluster( centers_fg, ...
                fcenters', length(fd), paramPropagate.nclus );
        end
    end
    centers_fg = centers_fg';
    % get the labels of boundary superpixels in each cluster
    bound_clus=cell(1,paramPropagate.nclus);   
    for i=1:length(fd)    
      bound_clus{featLabel_fg(i)}=[bound_clus{featLabel_fg(i)},fd(i)];  
    end
    for i=1:clus_num
        fcenters = centers_fg(i,:);
        f_x = repmat(fcenters(:,1)',[size(bound_clus{i},1) 1]);
        f_y = repmat(fcenters(:,2)',[size(bound_clus{i},1) 1]); 
        f_disxs = (superpixels.centres(bound_clus{i}(:),1) -f_x).^2;%%%%%%%%%%%%xÊÇ·´µÄ
        f_disys = (superpixels.centres(bound_clus{i}(:),2) -f_y).^2;
        f_disx =  min((3*mean(f_disxs)).^0.5,0.3*height);
        f_disy =  min((3*mean(f_disys)).^0.5,0.3*width);
        if f_disx<f_disy
            f_disy =  min(f_disy,f_disx*3);
        else
            f_disx =  min(f_disx,f_disy*3);
        end
        leftx = max(round(fcenters(:,1)-f_disx),1);
        rightx = min(round(fcenters(:,1)+f_disx),height);
        lefty = max(round(fcenters(:,2)-f_disy),1);
        righty = min(round(fcenters(:,2)+f_disy),width);
        mask(leftx:rightx,lefty:righty)=1;
    end
    mask(1:10,:)=0;mask(end-10:end,:)=0;mask(:,1:10)=0;mask(:,end-10:end)=0;
    mask = logical(mask);