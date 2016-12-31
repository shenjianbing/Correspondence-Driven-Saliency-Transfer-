%%%%%%%%%%%%%%计算mask%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [BSsaliency mask leftx rightx lefty righty]= centerBias(IniSal_region,IniSal,superpixels) 
[height,width] = size(IniSal);
nLabel = double(max(superpixels.Label(:)));
Label = superpixels.Label;
L{1} = uint32(Label);
bdIds=[Label(1,:)';Label(end,:)';Label(:,1);Label(:,end)];
bdIds = unique(bdIds);

colDistM = squareform(pdist(superpixels.Lab));

[conSPix Conedge]= find_connect_superpixel2( Label, nLabel, height ,width );  
ConSPix=sparse([Conedge(:,1);Conedge(:,2)],[Conedge(:,2);Conedge(:,1)], ...
         [ones(size(Conedge(:,1)));ones(size(Conedge(:,1)))],nLabel,nLabel);
ConSPix = full(ConSPix);
ConSPix = ConSPix -diag(diag(ConSPix));
clipVal = EstimateDynamicParas(ConSPix, colDistM);

ConSPix = ConSPix+eye(size(ConSPix));
              
b_geoDist = GeodesicSaliency(ConSPix, double(bdIds), colDistM, clipVal,IniSal_region(bdIds));
R = b_geoDist'.*IniSal_region;

R = (R-min(R))/(max(R)-min(R));
IniSal = reshape(superpixel2pixel(double(superpixels.Label),R),height,width);
IniSal = IniSal.*saliencycompact(IniSal);
S{1} = repmat(IniSal(:),[1 3]);
[ R, ~, ~ ] = getSuperpixelStats(S(1:1),L, double(max(Label(:))));
R = double(R(:,1));
R = (R-min(R))/(max(R)-min(R));
IniSal = reshape(superpixel2pixel(double(Label),R),height,width);

mask = (IniSal>mean(IniSal(:)));
mask = int32(mask>0);

fd = mask.*Label;
fd = unique(fd(:));
fd(fd==0) = [];
fd = sort(fd);
bd = int32(IniSal<=mean(IniSal(:))).*Label;
bd = unique(bd(:));
bd(bd==0) = [];

% colDistM = squareform(pdist(superpixels.Lab));
% colDistMx = colDistM;
% colDistMx(bd,:) = colDistMx(bd,:);
% colDistMx(:,bd) = colDistMx(:,bd);
% [conSPix Conedge]= find_connect_superpixel( Label, nLabel, height ,width );  
% ConSPix=sparse([Conedge(:,1);Conedge(:,2)],[Conedge(:,2);Conedge(:,1)], ...
%          [ones(size(Conedge(:,1)));ones(size(Conedge(:,1)))],nLabel,nLabel);
% ConSPix = full(ConSPix);
% ConSPix = ConSPix +eye(size(ConSPix));
       
        
% bcenters = superpixels.centres( bd,:);
% b_x = repmat(bcenters(:,1)',[nLabel 1]);
% b_y = repmat(bcenters(:,2)',[nLabel 1]);        
% x = repmat(superpixels.centres(:,1),[1 size(bd,1)]);
% y = repmat(superpixels.centres(:,2),[1 size(bd,1)]);        
% b_dis = ((x-b_x).^2+(y-b_y).^2).^0.5;
% b_dis = min(b_dis');
% b_dis = b_dis(:,fd);
% b_dis = (b_dis-min(b_dis(:)))/(max(b_dis(:))-min(b_dis(:)));
% b_dis = exp(2*b_dis);

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
        f_disxs = (superpixels.centres(bound_clus{i}(:),1) -f_x).^2;%%%%%%%%%%%%x是反的
        f_disys = (superpixels.centres(bound_clus{i}(:),2) -f_y).^2;
        f_disx =  (3*mean(f_disxs)).^0.5;
        f_disy =  (3*mean(f_disys)).^0.5;
        leftx = max(round(fcenters(:,1)-f_disx),1);
        rightx = min(round(fcenters(:,1)+f_disx),height);
        lefty = max(round(fcenters(:,2)-f_disy),1);
        righty = min(round(fcenters(:,2)+f_disy),width);
        mask(leftx:rightx,lefty:righty)=1;

    end
    mask(1:10,:)=0;mask(end-10:end,:)=0;mask(:,1:10)=0;mask(:,end-10:end)=0;
    mask = logical(mask);
    
    BSsaliency = IniSal;
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% f_x = repmat(fcenters(:,1)',[nLabel 1]);
% f_y = repmat(fcenters(:,2)',[nLabel 1]);        
% x = repmat(superpixels.centres(:,1),[1 size(fd,1)]);
% y = repmat(superpixels.centres(:,2),[1 size(fd,1)]);        
% f_dis = ((x-f_x).^2+(y-f_y).^2).^0.5;
% f_dis = f_dis(fd,:);
% f_dis = sum(f_dis,2)/size(fd,1);
% f_dis = (f_dis-min(f_dis(:)))/(max(f_dis(:))-min(f_dis(:)));
% f_dis = exp(5*f_dis);
% 
% border =  unique([find(fcenters(:,1)<height/4);find(fcenters(:,1)>3*height/4);find(fcenters(:,2)<width/4);find(fcenters(:,2)>3*width/4)]);
%         
% b_geoDist = GeodesicSaliency(ConSPix, double(bd), colDistM, 0,false,[]);
% b_geoDist = b_geoDist / max(b_geoDist(:)); 
% b_geoDist(:,bd) = [];
% b_geoDist = exp(b_geoDist);
% adjcMatrix_virtual = tril(ConSPix, -1);
% edgeWeight = colDistMx(adjcMatrix_virtual > 0);
%     
% f_geoDist = [];
% for i = 1:size(fd,1)
%     f_geoDist = [f_geoDist;graphshortestpath(sparse(adjcMatrix_virtual), fd(i), 'directed', false, 'Weights', edgeWeight)];             
% end
% f_geoDist(:,bd) = [];
% %geoDistMatrix = graphallshortestpaths(sparse(adjcMatrix), 'directed', false, 'Weights', edgeWeight);
% measure = max(f_geoDist).*f_dis'./b_geoDist./b_dis;
% if size(border,1)>1
%     measure(border) = measure(border)*2;
% end
% [~,pos] = min(measure);
% %[~,pos] = max(dis);
% center = fd(pos);
% 
% mask = zeros(height,width);
% % mask(Label==center) = 1;
% 
% 
% fcenters = superpixels.centres( center,:);
% f_x = repmat(fcenters(:,1)',[size(fd,1) 1]);
% f_y = repmat(fcenters(:,2)',[size(fd,1) 1]); 
% f_disxs = (superpixels.centres(fd(:),1) -f_x).^2;%%%%%%%%%%%%x是反的
% f_disys = (superpixels.centres(fd(:),2) -f_y).^2;
% 
% ratio = 4:-0.5:2;
% 
% for i = ratio
% f_disx =  (i*mean(f_disxs)).^0.5;
% f_disy =  (i*mean(f_disys)).^0.5;
% leftx = max(fcenters(:,1)-f_disx,1);
% rightx = min(fcenters(:,1)+f_disx,height);
% lefty = max(fcenters(:,2)-f_disy,1);
% righty = min(fcenters(:,2)+f_disy,width);
% if (rightx-leftx)*(righty-lefty)/(height*width)<0.8
%     break;
% end
% end
% mask(:) = 0;
% mask(leftx:rightx,lefty:righty)=1;
% mask = logical(mask);
% 
% [X,Y] = meshgrid(1:width,1:height); 
% BSsaliency = (exp(-((Y-fcenters(:,1)).^2/(f_disx^2)+ (X-fcenters(:,2)).^2/(f_disy^2))));
% BSsaliency = double(BSsaliency/max(BSsaliency(:)));





