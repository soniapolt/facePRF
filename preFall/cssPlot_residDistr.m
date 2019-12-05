% loads & plots SEM of beta estimates across ROIS

clear all; close all;

subjs = {'TH'};%
task = '';
expt = 'fixPRF';

saveFig = 1;

minR2 = 20;          % cutoff for vox selection
ROIs= standardROIs;

whichStim = 'photo';%'internal';%
whichModel = 'kayCSS';%'cssExpN';%'cssShift';%
whichM = 3; % 1 = mean, 2 = mode/peak, 3 = median
hems = {'rh' 'lh'};
fitSuffix = '';%'_orig';%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 11; titleSize = 14;
nBins = 20; % histogram bins

% now we load in the data from both hemispheres, and threshold across
load(pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems,fitSuffix,task));
ROInum = cellNum(ROIs,info.ROIs);
subjNum = cellNum(subjs,info.subjs);

if length(subjNum) == 1 roi = subj(subjNum).roi; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure: SEMs of betas for this ROI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

titleText = [whichModel ' beta SEMs, Subj: '];
titleText = [titleText strTogether(subjs) ' (voxels R^2 > ' num2str(minR2) '), ' whichStim ' stim, ' whichModel ' model'];
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if saveFig && onLaptop figSize = [0 0 1 1]; else figSize = [.2 .1 8 .8]; end
    niceFig(figSize,fontSize);
    numPlots = [2 ceil(length(ROIs)/2)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for r = 1:length(ROIs)
    subplot(numPlots(1),numPlots(2),r)
    
    for c = 1:length(roi(1).fits)
        cSEMs{c} = vertcat(roi(ROInum(r)).fits(c).vox.sems)'; 
        %cSEMs{c} = mean([roi(ROInum(r)).fits(c).vox.sems]); 
    end
    
%     clear h; clear p; clear stats;
%     for v = 1:length(roi(ROInum(r)).fits(1).vox)
%      [h(v),p(v),~,stats(v)] = ttest2([roi(ROInum(r)).fits(1).vox(v).sems],[roi(ROInum(r)).fits(2).vox(v).sems]);
%     end
    
    %niceHist(p,condColors(4),1);
    
    plotDistr(cSEMs,1,{roi(ROInum(r)).fits.cond},nBins,whichM,1);
    xlabel([ROIs{r} ': '  num2str(length(roi(ROInum(r)).fits(1).vox)) ' voxels'],'fontSize',titleSize,'interpreter','none','FontWeight','bold');
%    xlabel([ROIs{r} ': ' num2str(sum(h)) '/' num2str(v) ' signif. voxels'],'fontSize',titleSize,'interpreter','none','FontWeight','bold');
    
    end
    superTitle(titleText,titleSize,.025);

if saveFig == 1
    if length(subjs) == 1 txt = ['subj' subjs{1}]; else txt = ['groupN' num2str(length(subjs))]; end
    if length(hems) == 1
        txt = [txt '_' hems{1}]; end
        txt = [whichModel '_' txt];
    niceSave([dirOf(pwd) 'figures/' expt '/sem/'],txt); % just save pngs, since these can be generated pretty quickly
end

if onLaptop playSound; end