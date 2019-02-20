% loads & plots distributions of XY changes/anything else per each subject

clear all; close all;

subjs = {'SP' 'TH' 'DF' 'EM' 'MG' 'JG'};% 
%task = 'fix';
expt = 'fixPRF';

saveFig = 1;

minR2 = 20;          % cutoff for vox selection
ROIs= standardROIs;%{'IOG_faces' 'pFus_faces' 'mFus_faces'};%'V1' 'V2' 'V3' 'hV4'

whichStim = 'photo';%'internal';
whichModel = 'kayCSS';%'cssExpN';%'cssShift';%
whichM = 3; % 1 = mean, 2 = mode/peak, 3 = median
hems = {'rh' 'lh'};
fitSuffix = '';%'_orig';%

distLong = {'Eccen (dva)' 'Size (2*SD/sqrt(N)) (dva)' 'R2'};
distShort = {'eccen' 'size' 'r2'};

switch expt
    case 'fixPRF' 
        baseCond = [2]; compConds = [1]; 
    case 'compPRF' 
        baseCond = [3]; compConds = [1 2]; 
end
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
% figure 1: scatterplot of parameters for this ROI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p = 1:length(distShort)
    
    titleText = [expt ' ' distLong{p} ', Subj: '];
    titleText = [titleText strTogether(subjs) ' (voxels R^2 > ' num2str(minR2) '), ' whichStim ' stim, ' whichModel ' model'];
    for c = 1:length(compConds)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if saveFig && onLaptop figSize = [0 0 1 1]; else figSize = [.2 .1 .8 .8]; end
    niceFig(figSize,fontSize);
    numPlots = [2 ceil(length(ROIs)/2)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for r = 1:length(ROIs)
        subplot(numPlots(1),numPlots(2),r)
        if length(subjNum) == 1 prfs = subj(subjNum).roi(ROInum(r)); else prfs = roi(ROInum(r)); end
        eval(['bPars = [prfs.fits(baseCond).vox.' distShort{p} '];']);
        eval(['cPars = [prfs.fits(compConds(c)).vox.' distShort{p} '];']);
        
        scatterCent(bPars,cPars,condColors(compConds(c)),...
                prfs.fits(baseCond).cond,prfs.fits(compConds(c)).cond,distLong{p},fontSize);

        title([ROIs{r} ' (' num2str(length(roi(r).fits(1).vox)) ' vox)'],'fontSize',titleSize,'interpreter','none','FontWeight','bold');
        
    end
    superTitle(titleText,titleSize,.97);
    
    if saveFig == 1
        if length(subjs) == 1 txt = ['scatter_v' prfs.fits(c).cond '_subj' subjs{1}]; else txt = ['groupN' num2str(length(subjs))]; end
        if exist('task','var') [distShort{p} '_task' task '_' txt ];
        else txt = [distShort{p} '_' txt ];
        end
        if length(hems) == 1
            txt = [txt '_' hems{1}]; end
        
        if ~containsTxt(whichStim,'photo')
            txt = [whichStim '_' txt];  end
        
        txt = [whichModel '_' txt];
        
        niceSave([dirOf(pwd) 'figures/' expt '/params/'],txt); % just save pngs, since these can be generated pretty quickly
    end
    end
end
if onLaptop playSound; end