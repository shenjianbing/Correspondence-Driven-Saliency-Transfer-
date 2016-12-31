clc
clear
addpath( genpath( '.' ) );

%database = 'MSRA10K';
%database = 'ECSSD';
database = 'DUT';
%database = 'MSRA5000';
%database = 'PASCAL';
setting = setEnvironment(database);
optPara.range =1:setting.para.nTest;
load (fullfile(setting.path.database, 'splitTrainingTest.mat'));
for m = optPara.range
    [im1 imgName1]=STimread(setting.path.annotation,testndx(m));
    %[im1 imgName1]= STimread(setting.path.imgs,testndx(m));    
     imwrite(im1(:,:,1),[imgName1 '.png']);
end
