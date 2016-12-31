clc
clear
addpath( genpath( '.' ) );

%database = 'MSRA10K';
database = 'ECSSD';
%database = 'MSRA1000';
%database = 'DUT';
setting = setEnvironment(database);

optPara.nCandidates = setting.nCandidates;
optPara.prior_weight =0.06;
optPara.spatial_weight =20;
optPara.isDisplay =true;
optPara.isPrint =true;
optPara.isSaveEps =true;

optPara.nNeighbors =setting.K;
optPara.alpha =setting.SIFTflowpara.alpha;
optPara.KNNmethod =setting.KNNmethod;
optPara.epsilon =10000;
optPara.range =1:setting.para.nTest;
optPara.trainingRatio =1;


%load data
load (fullfile(setting.path.database, 'splitTrainingTest.mat'));
gist = STgist(setting.path.imgs, setting.GISTparam, setting.path.gist);
gistTest = gist(testndx, :);
gistTraining = gist(trainingndx, :);
MAEmin = 10;
    for m = optPara.range(1:end)
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
    K = max(setting.nCandidates,sum(dist<=dist(1)*optPara.epsilon));
    knn = knn(1:K);
    dist = dist(1:K);
%     figure; 
%     imshow(STimread(setting.path.imgs,testndx(m)));
%     for i = 1:K
%         figure;
%         imshow(STimread(setting.path.imgs,trainingndx(knn(i))));
%     end
    [im1, imgName1]= STimread(setting.path.imgs,testndx(m));
    [height,width,~] = size(im1);
     disp(imgName1); 
    load(fullfile(setting.path.matching_saliency.home,[imgName1 '.mat']));
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
    UMsaliency = UMsaliency/max(UMsaliency(:));    
    imwrite(UMsaliency,[setting.path.saliency '\' imgName1 '_bp.png']);
   
    end

