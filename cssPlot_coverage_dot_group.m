% loads & plots some basics re: invPRF fits (css model)
% SP 10/25/17
% _indiv code allows us to plot individual ROIs in multiple plots
% 10/8/18 edit: made more flexible across expriment with different numbers
% of fit conditions

clear all; close all;

subjs = prfSubjs;%{'george'};%{'SP' 'DF' 'EM' 'TH' 'MG' 'JG'};%;{'DF'};%
expt = 'fixPRF';%'nhp';%
plotRect = 3.2;%4.0719;    % plotCoverage has support for plotCirc; this adds monkey imsize to our data

minR2 =20;          % cutoff for vox selection
ROIs= {'V1'};%['hV4' standardROIs('face')];% %{'PL' 'ML'};%
sampleVox = 300; % how many randomly selected voxels are we plotting?
plotSize = 1;

saveFig = 1;

whichStim = 'outline';%'photo';%'binary';%'eyes';%'internal';%
whichModel = 'kayCSS';%''cssExpN';%cssShift';%
fitSuffix = '';

hems = {'rh' 'lh'};%{''};% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 12; titleSize = 14;

% now we load in the data from both hemispheres, and threshold across
load(pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems,fitSuffix));
ROInum = cellNum(ROIs,info.ROIs);
subjNum = cellNum(subjs,info.subjs);

for r = 1:length(ROIs)
    if length(subjNum) == 1 fits = subj(subjNum).roi(ROInum(r)).fits;
    else fits =roi(ROInum(r)).fits; end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  create supertitle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    titleText = [whichModel ' ' hemText(hems) ' ' ROIs{r} ', Subjs: ' strTogether(subjs) ...
        ' (' num2str(length(fits(1).vox)) ' voxels R^2 > ' num2str(minR2) ', sampleVox = ' num2str(sampleVox) '), ' whichStim];
    
    if onLaptop niceFig([.1 .1 .8 .8],fontSize); else
        niceFig([0 .2 .6 .25*(3)],fontSize); end
    numPlots = [1 length(fits)]; pl = 1;
   
    %input('Press A Key');
    pause(3);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1) coverage, 
    for c = 1:length(fits)
        %niceFig([.1 .1 1 1],16)
    subplot(numPlots(1),numPlots(2),pl)
    % plotCoverage(vox,color,leg,ppd,res,plotSize,alphaGain,sampleVox,centerMass,plotCirc)
    plotCoverage(fits(c).vox,condColors(r),'',roi(1).fits(1).ppd,roi(1).fits(1).res,plotSize,1,sampleVox);%,1,plotRect/2);
    
%     if strcmp(expt,'nhp') && plotRect>0
%     hold on; rectangle('Position',[-plotRect/2,-plotRect/2,plotRect,plotRect]); end

    t = title(roi(1).fits(c).cond);
    set(t,'visible','on');
    pl = pl+1;
    end
    
    subplotresize(numPlots(1), numPlots(2),.7,.7);
    
    superTitle(titleText,titleSize,.05);
    if saveFig
        txt = [hemText(hems) '_' ROIs{r} '_' whichModel  '_' whichStim ];
        if plotSize txt = [txt '_sizeLines']; end
        if ~isequal(minR2,20)
            txt = [txt '_r2-' num2str(minR2)];end
        niceSave([dirOf(pwd) 'figures/' expt '/lineCoverage/'],txt,[],subjs); % just save pngs, since these can be generated pretty quickly
    end
end % ROIs

if onLaptop playSound; end