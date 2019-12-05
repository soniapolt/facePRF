% s1_calcAveragedCoverage.m
%
% This script will loop through subject retinotopy data and
% will then create averaged coverage maps
%
% Adapted from DF 10/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

subjs = prfSubjs;%{'TH' 'DF' 'EM' 'JG' 'MG' 'SP'};
expt = 'fixPRF';

minR2 = 20;          % cutoff for vox selection
ROIs= standardROIs('face');
fitSuffix = '';

whichStim = 'photo';
whichModel = 'kayCSS';
hems = {'lh' 'rh'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems,fitSuffix));

ROInum = cellNum(ROIs,info.ROIs);
subjNum = cellNum(subjs,info.subjs);

exptDir = fullfile(raid,expt);
savePath = [exptDir '/coverage']; checkDir(savePath);

pl = 1;
for c = 1:2
for r = 1:length(ROInum)
   subplot(2,length(ROIs),pl)
   mrvCoverage(roi(ROInum(r)).fits(c).vox,[hemText(hems) '_' ROIs{r}]);
end
end

    %Save the coverage info 
    % saveFile = fullfile(savePath, [hems{h} '_coverage_data_ve', num2str(ve_cutoff*100)]);
   