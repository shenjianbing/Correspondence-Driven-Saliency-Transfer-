clc
clear
addpath( genpath( '.' ) );
database = 'subDUT';
fprintf('Creating dataset.\n');
setting=setEnvironment(database);

Imgs = imdir(setting.path.imgs);
Nimages = length(Imgs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute GIST
gist = STgist(setting.path.imgs, setting.GISTparam, setting.path.gist);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Split dataset in test and training and store the indices in the root
% folder for the dataset

testndx = fix(linspace(1, Nimages, setting.para.nTest));
trainingndx = setdiff(1:Nimages, testndx);
save(fullfile(setting.path.database, 'splitTrainingTest.mat'), 'testndx', 'trainingndx') 

fprintf('Dataset has been created.\n');



    

