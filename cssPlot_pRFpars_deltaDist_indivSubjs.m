% loads & plots distributions of XY changes/anything else per each subject

clear all; close all;

subjs = prfSubjs;
expt = 'fixPRF';

saveFig = 0;
convertDVA = 1;

minR2 = 'r2-20';          % cutoff for vox selection
ROI = standardROIs(7);%['V1' standardROIs('face')]

whichStim = 'outline';%'photo';%'internal';%'eyes';%
whichModel = 'kayCSS';%'inflipCSSn';%'flipCSSn';%'cssExpN';%'cssShift';%
whichM = 3; % 1 = mean, 2 = mode/peak, 3 = median
mS = {'mean' 'mode' 'median'};
plotPars = {'Y'};%{'gain' 'r2' 'Y' 'X' 'size'};%{'Ydeg'}%
parTitles = {'Y estim'};%{'Gain Estim' 'Estimated R^{2}' 'Y Estim.' 'X Estim.' 'Size [2xSD/sqrt(N)] (dva)'};

hems = {'lh' 'rh'};
fitSuffix = '';%'_orig';%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 11; titleSize = 14;
nBins = 20; % histogram bins

% now we load in the data from both hemispheres, and threshold across
load(pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems,fitSuffix));
ROInum = cellNum(ROI,info.ROIs);
subjNum = cellNum(subjs,info.subjs);

if length(subjNum) == 1 roi = subj(subjNum).roi; else subj = subj(subjNum); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1: distribution of parameters for this ROI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p = 1:length(plotPars)
    
    titleText = [whichModel ' ' parTitles{p} ', Subj: '];
    titleText = [titleText strTogether(subjs) ' (voxels R^2 > ' num2str(minR2) '), ' whichStim ' stim, ' whichModel ' model'];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if saveFig && onLaptop figSize = [.1 .1 .8 .8]; else figSize = [.2 .1 8 .8]; end
    f1 = niceFig(figSize,fontSize);
    numPlots = [2 ceil(length(subjs)/2)];pl = 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    sPars = []; mPars = [];
    for s = 1:length(subjs)
        
        subplot(numPlots(1),numPlots(2),pl)
        
        for c = 1:length(roi(1).fits)
            
            fits = subj(subjNum(s)).roi(ROInum).fits(c);
            
            try
                cPars{c} = getPar(plotPars{p},fits,1);
            catch cPars{c} = NaN; end % missing ROIs
            
            % aggregate across subjects
            sPars{s,c} = cPars{c}; % full distribution for this subject
            eval(['mPars(s,c) = nan' mS{whichM} '(cPars{c});']); % mean value for this subject
        end
        hues = [1 .5];
        if containsTxt(plotPars{p},'size')
            xl = [0 5] ;
        else xl = [-2 2]; end
        dPars{s} = sPars{s,1}-sPars{s,2};
        
        plotDistr({dPars{s}},1,{''},20);
        xlabel('Inverted minus Upright');
        set(gca,'TickDir','out');
        pl = pl+1;
        axis square;
        
        title({subjs{s};parTitles{p}});
    end
    
    superTitle(titleText,titleSize,.025);
    
    f2=figure;
    [h,meds, normcounts] = plotMeanDistr(dPars,nBins,roiColors(ROI),1); %xlim(xl);
    if saveFig == 1
        figure(f(1));
        txt = [plotPars{p} '_' ROI 'indivSubjs_' hemText(hems)];
        niceSave([dirOf(pwd) 'figures/' expt '/deltaDist/'],txt,[],[],{'png' 'svg'}); % just save pngs, since these can be generated pretty quickly
    end
end
if onLaptop playSound; end