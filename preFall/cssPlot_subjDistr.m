% loads & plots distributions of XY changes/anything else per each subject

clear all; close all;

subjs = {'MG' 'JG' 'TH' 'EM' 'DF' 'SP' 'JP' 'MH' 'JW' 'MN' 'JJ'};
task = '';
expt = 'fixPRF';

saveFig = 1;

minR2 = 20;          % cutoff for vox selection
ROIs= standardROIs('face');%{'IOG_faces' 'pFus_faces' 'mFus_faces'};%'V1' 'V2' 'V3' 'hV4' 
comps = [1 2];%; 1 3; 2 3];
whichPlot = 1; % 1 = subjs, 2 = group, 3 = both
fitSuffix = '';

whichStim = 'photo';
whichModel = 'kayCSS';%'cssShift';%

hems = {'rh' 'lh'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 11; titleSize = 14;
nBins = 20; % histogram bins
% now we load in the data from both hemispheres, and threshold across
load(pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems,fitSuffix,task));
ROInum = cellNum(ROIs,info.ROIs);
subjNum = cellNum(subjs,info.subjs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  create supertitle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

titleText = [expt ' (' hemText(hems) ', Subj: '];
titleText = [titleText strTogether(subjs) ' (voxels R^2 > ' num2str(minR2) '), ' whichStim ' stim, ' whichModel ' model'];

for cc = 1:size(comps,1)
    baseCond = comps(cc,1);
    c = comps(cc,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveFig == 1 figSize = [0 0 1 1]; else figSize = [.2 .1 1 .8]; end
niceFig(figSize,fontSize);
numPlots = [2 length(ROIs)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1: (repeat for each comparison)
% top row) XY shift magnitude along visual heirarchy
% bottom row) size change along visual heirarchy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for r = 1:length(ROIs)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1) shifts across ROIs
    subplot(numPlots(1),numPlots(2),r)
    for s = 1:length(subjs)
        % get XY shifts
        [~, subjCh,subjSign] = plotXYshift2(subj(subjNum(s)).roi(ROInum(r)).fits(baseCond).vox,subj(subjNum(s)).roi(ROInum(r)).fits(c).vox,roi(1).fits(1).ppd,roi(1).fits(1).res,1);
        % plot them
        subjData{s} = subjCh/roi(1).fits(1).ppd.*subjSign;
    end
    
    plotDistr(subjData,whichPlot,subjs,nBins);
    
    title(ROIs{r},'fontSize',titleSize,'interpreter','none','FontWeight','bold'); 
    xlabel(['(-) Towards Center, (+) Away From Center'],'FontSize',fontSize);
    
    if r == 1
    ylabel({[subj(s).roi(r).fits(baseCond).cond ' to ' subj(s).roi(r).fits(c).cond ];'XY Shifts (dva)'},'FontSize',titleSize+8,'FontWeight','bold');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2) size change across ROIs
    subplot(numPlots(1),numPlots(2),r+numPlots(2))
    for s = 1:length(subjs)
        % get XY shifts
        [sizeCh] = plotSizeChange(subj(subjNum(s)).roi(ROInum(r)).fits(baseCond).vox,subj(subjNum(s)).roi(ROInum(r)).fits(c).vox,roi(1).fits(1).ppd,roi(1).fits(1).res,1);
        % plot them
        subjData{s} = sizeCh/roi(1).fits(1).ppd;
    end
    
    plotDistr(subjData,whichPlot,subjs,nBins); 
    
    title(ROIs{r},'fontSize',titleSize,'interpreter','none','FontWeight','bold'); 
    xlabel(['(-) Smaller, (+) Bigger'],'FontSize',fontSize);
    
    if r == 1
    ylabel({[subj(s).roi(r).fits(baseCond).cond ' to ' subj(s).roi(r).fits(c).cond ];'Size Changes (dva)'},'FontSize',titleSize+8,'FontWeight','bold');
    end
    
    end
    superTitle(titleText,titleSize,.97);


if saveFig == 1
    if whichPlot == 1 txt = 'subjs'; elseif whichPlot == 2 txt = 'group'; else txt = []; end
    switch expt
            case 'invPRF3'
                txt = ['dist_' task 'Task_' fits(c).cond '_' txt '_'];
            case 'fixPRF'
                txt = ['dist_' roi(1).fits(baseCond).cond 'vs' roi(1).fits(c).cond '_' txt '_'];
            case 'compPRF'
                txt = ['dist_' roi(1).fits(baseCond).cond 'vs' roi(1).fits(c).cond '_' txt '_'];
    end
    if length(hems) == 1
        txt = [txt hems{1} '_']; end
    if ~containsTxt(whichStim,'photo')
            txt = [whichStim '_' txt];  end
    txt = [whichModel txt];
    niceSave([dirOf(pwd) 'figures/' expt '/subjDists/'],txt,ROIs); % just save pngs, since these can be generated pretty quickly
end

end