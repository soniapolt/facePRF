% loads & plots distributions of XY changes/anything else per each subject

clear all; close all;

subjs = {'TH' 'DF' 'EM' 'MG' 'JG' 'SP'};%
%task = 'fix';
expt = 'fixPRF';

saveFig = 0;

minR2 = 20;          % cutoff for vox selection
ROIs= {'hV4' 'IOG_faces' 'pFus_faces' 'mFus_faces'};%'V1' 'V2' 'V3' 'hV4' 

whichStim = 'internal';%'eyes';%'photo';%
whichModel = 'kayCSS';%'cssExpN';%'cssShift';%
whichM = 3; % 1 = mean, 2 = mode/peak, 3 = median
hems = {'rh' 'lh'};
fitSuffix = '';%'_orig';%

distLong = {'Eccen (dva)' 'Size (2*SD/sqrt(N)) (dva)' 'R2'};
distShort = {'eccen' 'size' 'r2'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 11; titleSize = 14;
nBins = 20; % histogram bins

% now we load in the data from both hemispheres, and threshold across
load(pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems,fitSuffix));
ROInum = cellNum(ROIs,info.ROIs);
subjNum = cellNum(subjs,info.subjs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1: distribution of parameters for this ROI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p = 1:length(distShort)
    
titleText = [expt ' ' distLong{p} ', Subj: '];
titleText = [titleText strTogether(subjs) ' (voxels R^2 > ' num2str(minR2) '), ' whichStim ' stim, ' whichModel ' model'];
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if saveFig && onLaptop figSize = [0 0 1 1]; else figSize = [.2 .1 8 .8]; end
    niceFig(figSize,fontSize);
    numPlots = [2 ceil(length(ROIs)/2)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for r = 1:length(ROIs)
    subplot(numPlots(1),numPlots(2),r)
    if length(subjNum) == 1 prfs = subj(subjNum).roi(ROInum(r)); else prfs = roi(ROInum(r)); end
    for c = 1:length(roi(1).fits)
        eval(['cPars{c} = [prfs.fits(c).vox.' distShort{p} '];']);
    end
    
    plotDistr(cPars,1,{prfs.fits.cond},nBins,whichM);
    
    title(ROIs{r},'fontSize',titleSize,'interpreter','none','FontWeight','bold'); 
    xlabel(distLong{p},'fontSize',titleSize);
    
    end
    superTitle(titleText,titleSize,.97);

if saveFig == 1
    if length(subjs) == 1 txt = ['subj' subjs{1}]; else txt = ['groupN' num2str(length(subjs))]; end
    switch expt
            case 'invPRF3'
                txt = [ distShort{p} '_' task 'Task_' txt ];
            case 'fixPRF'
                txt = [distShort{p} '_' txt ];
            case 'compPRF'
                txt = [distShort{p} '_' txt ];
    end
    if length(hems) == 1
        txt = [txt '_' hems{1}]; end
    %if ~containsTxt(whichStim,'photo')
        txt = [whichStim '_' txt]; %end
    
        txt = ['distr_' whichModel '_' txt];
    niceSave([dirOf(pwd) 'figures/' expt '/params/'],txt); % just save pngs, since these can be generated pretty quickly
end
end
if onLaptop playSound; end