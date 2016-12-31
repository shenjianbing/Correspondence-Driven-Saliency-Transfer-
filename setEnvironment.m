% procedure to set the environmental parameters for scene matching

function setting = setEnvironment(database)


setting.path.home  = pwd;
setting.database = database;
setting.path.database     = fullfile(setting.path.home,'Datasets',database);
setting.path.imgs       = fullfile(setting.path.database,'Imgs');
setting.path.annotation   = fullfile(setting.path.database,'Annotations');
setting.path.gist         = fullfile(setting.path.database,'Gist');
setting.path.saliency         = fullfile(setting.path.database,'saliency');
setting.path.matching.home = fullfile(setting.path.database,'Matching');
setting.path.matching_saliency.home = fullfile(setting.path.database,'Matching_saliency');
% recursively create the directories of the fields in setting.path
structMkdir(setting.path);
    
%--------------------------------------------------------------------
% set the number of test images
%--------------------------------------------------------------------
setting.para.nTest = 10;

%-------------------------------------------------------------
% SOME PARAMETERS
%-------------------------------------------------------------

% gist
setting.GISTparam.imageSize                    = 256;
setting.GISTparam.orientationsPerScale = [8 8 8 8];
setting.GISTparam.numberBlocks            = 4;
setting.GISTparam.fc_prefilt                      = 4;
setting.GISTparam.boundaryExtension = 32; % number of pixels to pad

% hog
setting.VWparamhog.imageSize               = [256 256]; % it works also with non-square images
setting.VWparamhog.grid_spacing         = 1; % distance between grid centers
setting.VWparamhog.patch_size              = 16; % size of patch from which to compute SIFT descriptor (it has to be a factor of 4)
setting.VWparamhog.NumVisualWords = 200; % number of visual words
setting.VWparamhog.Mw                            = 2; % number of spatial scales for spatial pyramid histogram
setting.VWparamhog.descriptor               = 'hog'; 
setting.VWparamhog.w                                = floor(setting.VWparamhog.patch_size/2*2.5); % boundary for HOG

% dense sift
setting.SIFTparam.grid_spacing = 1; % distance between grid centers
setting.SIFTparam.cell_size    = 3; % size of cell from which to compute SIFT descriptor

% the number of nearest neighbors of scene retrieval
setting.K                           = 30;
setting.nCandidates           = 10;
setting.saliencylevel           = 3;
setting.patchratio               = 0.25;%4*4
setting.valScale = 60;
setting.alpha = 0.0004;
% setting the KNN methods; 
setting.KNNmethod  = 'gistL2';
% setting.KNNmethod = 'gistCorr';
% setting.KNNmethod = 'hog';
% setting.KNNmethod = 'gt';
% setting.KNNmethod = 'all';

% parameters for SIFT flow
setting.SIFTflowpara.alpha          = 0.7*255;  
setting.SIFTflowpara.d              = setting.SIFTflowpara.alpha*20;
setting.SIFTflowpara.gamma          = 0.0001;
setting.SIFTflowpara.nlevels        = 4;
setting.SIFTflowpara.topwsize       = 20;
setting.SIFTflowpara.wsize          = 5;
setting.SIFTflowpara.nIterations    = 50; % increasing the number of iterations leads to higher accuracy yet more tme
setting.SIFTflowpara.nTopIterations = 100;

% SIFTflowpara.topwsize=20;
% SIFTflowpara.wsize=3;
% SIFTflowpara.nIterations=30;
% SIFTflowpara.nTopIterations=60;
