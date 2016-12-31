clc
clear
addpath( genpath( '.' ) );

%database = 'MSRA10K';
%database = 'ECSSD';
%database = 'MSRA1000';
%database = 'DUT';
database = 'subDUT';
%database = 'MSRA5000';

setting = setEnvironment(database);

optPara.nCandidates = setting.nCandidates;
optPara.nNeighbors =setting.K;
optPara.prior_weight =0.06;
optPara.spatial_weight =20;
optPara.isDisplay =true;
optPara.isPrint =true;
optPara.isSaveEps =true;

optPara.alpha =setting.SIFTflowpara.alpha;
optPara.KNNmethod =setting.KNNmethod;
optPara.epsilon =10000;
optPara.range =1:setting.para.nTest;
optPara.trainingRatio =1;

RWRalpha = 0.0004;
pca_basis = [];
sift_size = 4;
%load data
load (fullfile(setting.path.database, 'splitTrainingTest.mat'));
gist = STgist(setting.path.imgs, setting.GISTparam, setting.path.gist);
gistTest = gist(testndx, :);
gistTraining = gist(trainingndx, :);

for m = optPara.range
    fprintf('Computing Saliency: %d -th test image: ',m);
    clear idx ik knn KNN energylist;
    switch optPara.KNNmethod
        case 'gistL2'
            [knn,dist] = LMgistquery(gistTest(m,:), gistTraining);
        case 'gistCorr'
            [knn,dist] = LMgistquery(gistTest(m,:), gistTraining,1);
        case 'hog'
            [knn,dist] = LMhistintersectionquery(hogTest(m,:), hogTraining);
        case'obj'            
            [knn,dist] = LMhistintersectionquery(objTest(m,:), objTraining);
        otherwise
    end
    dist = dist(1:optPara.nNeighbors);
    K = max(optPara.nCandidates,sum(dist<=dist(1)*optPara.epsilon));
    knn = knn(1:K);
    dist = dist(1:K);
    [im1, imgName1]= STimread(setting.path.imgs,testndx(m));
    [height,width,~] = size(im1);
    disp(imgName1); 

    %load(fullfile(setting.path.matching_saliency.home,['all' imgName1 '.mat']),'batch_siftenergy','siftenergy');
%     batch_siftenergy = batch_siftenergy(1:optPara.nNeighbors,:);
%     siftenergy = siftenergy(1:optPara.nNeighbors,:);
    if exist(fullfile(setting.path.matching_saliency.home,[imgName1 '_' num2digits(optPara.nNeighbors,2) '_' num2digits(optPara.nCandidates,2) '.mat']),'file')
        load(fullfile(setting.path.matching_saliency.home,[imgName1 '_' num2digits(optPara.nNeighbors,2) '_' num2digits(optPara.nCandidates,2) '.mat']));
    else    
        MATCHING = fullfile(setting.path.matching.home, ['alpha=' num2str(setting.SIFTflowpara.alpha/255)]);    
    %%%%%%%%saliency via local match%%%%%%%%%%%%%%%
        if exist(fullfile(MATCHING,['batchmatching_' num2str(0.25) '_' num2digits(m,4) '.mat']),'file')~=2
            continue;
        end
        matching = load(fullfile(MATCHING,['batchmatching_' num2str(setting.patchratio) '_' num2digits(m,4) '.mat']));
        BatchMatching = matching.batch_matching;

        batch_im1 = imresize(im1, setting.patchratio, 'bilinear'); 
        [h,w,~] = size(batch_im1);
        % compute the dense SIFT feature of the query image
        %batch_sift1 = mexDenseSIFT(batch_im1,setting.SIFTparam.cell_size,setting.SIFTparam.grid_spacing);
        batch_sift1 = ExtractSIFT(batch_im1, pca_basis, sift_size);
        batch_sift1 = reshape(batch_sift1,[h*w,128]);
        batch_siftenergy = zeros(K,h*w);
        batch_warpGT = zeros(K,h*w);
        for k = 1:K        
            [im2, ~]=STimread(setting.path.imgs,trainingndx(knn(k)));        
            batch_im2 = imresize(im2, setting.patchratio, 'bilinear'); 

            vx=double(BatchMatching(k).vx);
            vy=double(BatchMatching(k).vy);            
            batch_warpI2 = warpImage(double(batch_im2),vx,vy);
            %batch_warpsift2 = mexDenseSIFT(batch_warpI2,setting.SIFTparam.cell_size,setting.SIFTparam.grid_spacing);
            batch_warpsift2 = ExtractSIFT(batch_warpI2, pca_basis, sift_size);
            batch_warpsift2 = reshape(batch_warpsift2,[h*w,128]);
            siftdis =  sum((batch_sift1  - batch_warpsift2).^2,2);
            batch_siftenergy(k,:) = siftdis;

            [GT2, ~]=STimread(setting.path.annotation,trainingndx(knn(k)));
            GT2 = GT2(:,:,1); 
            batch_GT2 = imresize(GT2, setting.patchratio, 'bilinear');       
            batch_warpGT2 = warpImage(double(batch_GT2),vx,vy)/255;
            batch_warpGT2(batch_warpGT2<0)=0;
            batch_warpGT(k,:) = batch_warpGT2(:);
        end
        [~,batch_idx]=sort(batch_siftenergy);
        batch_idx = batch_idx(1:optPara.nCandidates,:);
        [XX,~]=meshgrid(1:h*w,1:optPara.nCandidates);
        batch_Candidates = full(sparse(batch_idx(:),XX(:),1,K,h*w));
        batch_saliency = batch_Candidates.*batch_warpGT;
        batch_saliency = sum(batch_saliency);
        LMsaliency_voting = imresize(reshape(batch_saliency,[h,w]),[height width],'nearest');       
        LMsaliency = (LMsaliency_voting-min(LMsaliency_voting(:)))/(max(LMsaliency_voting(:))-min(LMsaliency_voting(:)));   
    %%%%%%%%saliency via global match%%%%%%%%%%%%%%%
        if exist(fullfile(MATCHING,['pixelmatching' num2digits(m,4) '.mat']),'file')~=2
            continue;
        end
        matching = load(fullfile(MATCHING,['pixelmatching' num2digits(m,4) '.mat']));
        PixelMatching = matching.pixelmatching;
        %idx = matching.idx;
        %sift1 = mexDenseSIFT(im1,setting.SIFTparam.cell_size,setting.SIFTparam.grid_spacing);
        sift1 = ExtractSIFT(im1, pca_basis, sift_size);
        siftenergy = zeros(K,1);
        Msaliency_voting = zeros(height,width);
            for k = 1:K        
                [im2, ~] =STimread(setting.path.imgs,trainingndx(knn(k)));
                vx=double(PixelMatching(k).vx);
                vy=double(PixelMatching(k).vy);
    %             flow(:,:,1) = vx;
    %             flow(:,:,2) = vy;
    %             flowimg = flowToColor(flow);
                warpI2 = warpImage(double(im2),vx,vy);
                %warpsift2 = mexDenseSIFT(warpI2,setting.SIFTparam.cell_size,setting.SIFTparam.grid_spacing);
                warpsift2 = ExtractSIFT(warpI2, pca_basis, sift_size);
                siftdis=  (sift1  - warpsift2).^2;
                siftenergy(k) = sum(siftdis(:));
            end
        [~,idx]=sort(siftenergy);
        GMsaliency_voting = zeros(size(im1,1),size(im1,2));
        for k = 1:optPara.nCandidates
           [GT2, ~] =STimread(setting.path.annotation,trainingndx(knn(idx(k)))); 
           GT2 = double(GT2(:,:,1)); 
           vx=double(PixelMatching(idx(k)).vx);
           vy=double(PixelMatching(idx(k)).vy);
           warpGT2 = warpImage(GT2,vx,vy)/255;
           warpGT2(warpGT2<0)=0;
           GMsaliency_voting = GMsaliency_voting+warpGT2;
        end
        GMsaliency = (GMsaliency_voting-min(GMsaliency_voting(:)))/(max(GMsaliency_voting(:))-min(GMsaliency_voting(:)));   
        %imwrite(GMsaliency,[setting.path.saliency '/' imgName1 '_' num2digits(optPara.nNeighbors,2) '_' num2digits(optPara.nCandidates,2) '_GM' '.png']);
        %imwrite(LMsaliency,[setting.path.saliency '/' imgName1 '_' num2digits(optPara.nNeighbors,2) '_' num2digits(optPara.nCandidates,2) '_LM' '.png']);
        save(fullfile(setting.path.matching_saliency.home,[imgName1 '_' num2digits(optPara.nNeighbors,2) '_' num2digits(optPara.nCandidates,2) '.mat']),'LMsaliency_voting','LMsaliency','GMsaliency_voting','GMsaliency');
    end
    imwrite(GMsaliency,[setting.path.matching_saliency '/' imgName1 '_' num2digits(optPara.nNeighbors,2) '_' num2digits(optPara.nCandidates,2) '_GM' '.png']);
    imwrite(LMsaliency,[setting.path.matching_saliency '/' imgName1 '_' num2digits(optPara.nNeighbors,2) '_' num2digits(optPara.nCandidates,2) '_LM' '.png']);
    [im,w]=removeframe(im1);
    im = im1(w(3):w(4),w(5):w(6),:);
    LMsaliency_voting = LMsaliency_voting(w(3):w(4),w(5):w(6));
    LMsaliency = LMsaliency(w(3):w(4),w(5):w(6));
    GMsaliency_voting = GMsaliency_voting(w(3):w(4),w(5):w(6));
    GMsaliency = GMsaliency(w(3):w(4),w(5):w(6));
    [height,width,~] = size(im);
    Img{1} = double(im);   

    PixNum = height*width;        
    Attr=[ height ,width, 500, 20, PixNum ];
    imgVecR = reshape( Img{1}(:,:,1)', PixNum, 1);
    imgVecG = reshape( Img{1}(:,:,2)', PixNum, 1);
    imgVecB = reshape( Img{1}(:,:,3)', PixNum, 1); 

    [ Label, Sup1, Sup2, Sup3, RegionNum ] = SLIC( double(imgVecR), double(imgVecG), double(imgVecB), Attr );

    superpixels.Label = int32(reshape(Label+1,width,height)');
    superpixels.Lab = [Sup1 Sup2 Sup3] ;    
    L{1} = uint32(superpixels.Label);
    [ superpixels.colours, superpixels.centres, superpixels.t ] = getSuperpixelStats(Img(1:1), L, RegionNum );    
    UMsaliency_voting = min(GMsaliency_voting,LMsaliency_voting);
    UMsaliency = GMsaliency.*LMsaliency;
    R = getSuperpixelMeanValues(UMsaliency,superpixels.Label,double(max(superpixels.Label(:))));    
    UMsaliency_region = normalize(R);   
    UMsaliency = R(superpixels.Label);

   %% compute center bais     
   [BSsaliency, bias_mask]= centerBias(UMsaliency_region,UMsaliency,superpixels);

   [bp_R,fp_R] = computepro(BSsaliency,UMsaliency_voting,bias_mask,superpixels);  


    colDistM = squareform(pdist(superpixels.Lab));
    Label = [superpixels.Label;superpixels.Label(1,:)];    
    Label = [Label Label(:,1)];
    [conSPix, Conedge]= find_connect_superpixel2(Label, RegionNum, height+1 ,width+1 ); 

    valDistances=sqrt(sum((superpixels.Lab(Conedge(:,1),:)-superpixels.Lab(Conedge(:,2),:)).^2,2));  
    valDistances=normalize(valDistances);
    weights=exp(-setting.valScale*valDistances)+ 1e-5;
    weights=sparse([Conedge(:,1);Conedge(:,2)],[Conedge(:,2);Conedge(:,1)], ...
            [weights;weights],RegionNum,RegionNum);
    E = sparse(1:RegionNum,1:RegionNum,ones(RegionNum,1)); iD = sparse(1:RegionNum,1:RegionNum,1./sum(weights));
    P = iD*weights;

    RWR_R = (E-(1-RWRalpha)*P)\[bp_R fp_R]; %RWR
    RWR_R = RWR_R(:,2)./(RWR_R(:,1)+RWR_R(:,2));
    RWR_R = (RWR_R-min(RWR_R))/(max(RWR_R)-min(RWR_R));

    RWR=RWR_R(superpixels.Label);
    radius = floor(50*sqrt(mean(RWR(:))));
    RWR = morphSmooth(RWR, max(radius, 3));
    RWR = enhanceContrast(RWR, 10);
     
    RWRsaliency = zeros(w(1),w(2));
    RWRsaliency(w(3):w(4),w(5):w(6)) = RWR;
    UMMsaliency = zeros(w(1),w(2));
    UMMsaliency(w(3):w(4),w(5):w(6)) = UMsaliency;
    UMsaliency = UMMsaliency;
    BSSsaliency = zeros(w(1),w(2));
    BSSsaliency(w(3):w(4),w(5):w(6)) = BSsaliency;
    BSsaliency = BSSsaliency;
    
    imwrite(RWRsaliency,[setting.path.saliency '/' imgName1 '_' num2digits(optPara.nNeighbors,2) '_' num2digits(optPara.nCandidates,2) '.png']);
    imwrite(BSsaliency,[setting.path.saliency '/' imgName1 '_' num2digits(optPara.nNeighbors,2) '_' num2digits(optPara.nCandidates,2) '_BS' '.png']);
    imwrite(UMsaliency,[setting.path.saliency '/' imgName1 '_' num2digits(optPara.nNeighbors,2) '_' num2digits(optPara.nCandidates,2) '_UM' '.png']);
end




