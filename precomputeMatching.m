% set the environmental parameters
% load the training and test data
clc
clear
database = 'subDUT';
addpath( genpath( '.' ) );
setting = setEnvironment(database);

pca_basis = [];
sift_size = 4;
% range is used for parallelization 
range             = 1:setting.para.nTest;
K                 = setting.K ;
KNNmethod         = setting.KNNmethod ;

% if IsBackward==true, then we also estimate the flow from the each of the
% nearest neighbors back to the input image
IsBackward = false;
% if IsOverwrite==true, then the flow data will be overwritten anyway; otherwise, the computation will skip on existin
% files
IsOverwrite = true;

fprintf('Pre-computing bach-/pixel-level matching.\n');
%--------------------------------------------------------------------------------------------------------
%   create matching directory to save the flow data
%--------------------------------------------------------------------------------------------------------
MATCHINGPATH = fullfile(setting.path.matching.home,['alpha=' num2str(setting.SIFTflowpara.alpha/255 )]);
if exist(MATCHINGPATH,'dir')~=7
    mkdir(MATCHINGPATH);
end


%--------------------------------------------------------------------------------------------------------
% set multi thread
%--------------------------------------------------------------------------------------------------------
% if setting.IsMultiThread
%     nCores = setMultiThread;
% else
    nCores = 1;
%end
load (fullfile(setting.path.database, 'splitTrainingTest.mat'));
gist = STgist(setting.path.imgs, setting.GISTparam, setting.path.gist);
gistTest = gist(testndx, :);
gistTraining = gist(trainingndx, :);
%--------------------------------------------------------------------------------------------------------
% loop over the range
%--------------------------------------------------------------------------------------------------------
%
for m=range;
     fprintf('\tProcessing Patch-level Matching: %d -th test image',m);
%     if IsOverwrite == false && exist(fullfile(MATCHINGPATH,['matching' num2digits(m) '.mat']),'file') == 2
%         fprintf('SIFT matching for test image %d exists!\n',m);
%         continue;
%     end
%     fprintf('Find SIFT matching for test image %d ',m);
    
    %--------------------------------------------------------------------------------------------------------
    % retrieve nearest neighbors using a variety of methods
    %--------------------------------------------------------------------------------------------------------
    KNN = [];
    % use GIST L2 to retrieval K-nearest neighbors
    if strcmpi(KNNmethod,'gistL2') || strcmpi(KNNmethod,'all')
        knn = LMgistquery(gistTest(m,:), gistTraining,0);
        knn = knn(1:K);
        KNN = knn;
    end
    % use GIST normlized correlation to retrieval K-nearest neighbors
    if strcmpi(KNNmethod,'gistCorr') || strcmpi(KNNmethod,'all')
        knn = LMgistquery(gistTest(m,:), gistTraining,1);
        knn = knn(1:K);
        KNN = union(knn,KNN);
    end
    % use HOG histogram intersection to retrieve K-nearest neighbors
    if strcmpi(KNNmethod,'hog') || strcmpi(KNNmethod,'all')
        knn = LMhistintersectionquery(hogTest(m,:), hogTraining);
        knn = knn(1:K);
        KNN = union(knn,KNN);
    end    
    % use the ground-truth segmentation histogram intersection to retrieve K-NN
    if strcmpi(KNNmethod,'gt') || strcmpi(KNNmethod,'all')
        knn =  LMhistintersectionquery(objTest(m,:), objTraining);
        knn = knn(1:K);
        KNN = union(knn,KNN);
    end    
    clear knn;
    
    % load the image
    [im1 ~]= STimread(setting.path.imgs,testndx(m));
    batch_im1 = imresize(im1, 0.25, 'bilinear'); 
    % compute the dense SIFT feature of the query image
%     sift1 = mexDenseSIFT(im1,setting.SIFTparam.cell_size,setting.SIFTparam.grid_spacing);
%     batch_sift1 = mexDenseSIFT(batch_im1,setting.SIFTparam.cell_size,setting.SIFTparam.grid_spacing);
    sift1 = ExtractSIFT(im1, pca_basis, sift_size);
    batch_sift1 = ExtractSIFT(batch_im1, pca_basis, sift_size);
 % load the image
%     foo = Dtest(m).annotation;
%     im1 = imread(fullfile(fullfile(HOMEIMAGES,foo.folder),foo.filename));
%     im1 = cropImage(im1,SIFTparam.patch_size);
%     
%     subplottight(6,6,1);
%     imshow(showColorSIFT(sift1));
    fprintf('K=%d\n',length(KNN));
    BatchMatching = cell(1,length(KNN));
    
    for j = 1:length(KNN)
        [im2 ~]=STimread(setting.path.imgs,trainingndx(KNN(j)));
        batch_im2 = imresize(im2, 0.25, 'bilinear'); 
        % compute the dense SIFT feature of the nearest neighbor
        %batch_sift2 = mexDenseSIFT(batch_im2,setting.SIFTparam.cell_size,setting.SIFTparam.grid_spacing);
        batch_sift2 = ExtractSIFT(batch_im2, pca_basis, sift_size);
      
        % perform forward matching
        fprintf('\tMatching to No.%d NN (forward)...',j);

        [batch_vx,batch_vy] = DSPMatch(batch_sift1,batch_sift2); 
        BatchMatching{j}.vx = int16(batch_vx);
        BatchMatching{j}.vy = int16(batch_vy);
%         BatchMatching{j}.energylist = batch_energylist;
%         BatchMatching{j}.objmin = batch_energylist(1).data(end);
%         batch_energylist_all(j) = batch_energylist(1).data(end);
        batch_matching(j) = BatchMatching{j};
    end
    save(fullfile(MATCHINGPATH,['batchmatching_' num2str(0.25) '_' num2digits(m,4) '.mat']),'batch_matching','KNN');
    
    fprintf('\tProcessing Pixel-level Matching: %d -th test image',m);
    %PixelMatching = cell(1,setting.nCandidates);
    PixelMatching = cell(1,K);

    for j=1:K
        %index=KNN(idx(j));
        index=KNN(j);
        [im2 ~]=STimread(setting.path.imgs,trainingndx(index));       
        % compute the dense SIFT feature of the nearest neighbor
        %sift2 = mexDenseSIFT(im2,setting.SIFTparam.cell_size,setting.SIFTparam.grid_spacing);
        sift2 = ExtractSIFT(im2, pca_basis, sift_size);
        % perform forward matching
        fprintf('\tMatching to No.%d NN (forward)...',j);
        
        %[vx,vy,energylist]=SIFTflowc2f(sift1,sift2,setting.SIFTflowpara);
        [vx,vy] = DSPMatch(sift1,sift2); 
        PixelMatching{j}.vx = int16(vx);
        PixelMatching{j}.vy = int16(vy);
%         PixelMatching{j}.energylist = energylist;
%         PixelMatching{j}.objmin = energylist(1).data(end);

       % warpI2 = warpImage(im2double(im2),vx,vy)/255;
%     end
%     for j = 1:setting.nCandidates
        pixelmatching(j) = PixelMatching{j};
    end
    save(fullfile(MATCHINGPATH,['pixelmatching' num2digits(m,4) '.mat']),'pixelmatching');
    
end