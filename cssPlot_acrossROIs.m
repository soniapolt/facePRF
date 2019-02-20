% loads & plots changes along the visual heirarchy for cssFit properties
% currently - across all subject in pRFset file

clear all; close all;

task = '';
expt = 'fixPRF';%'invPRF3';%

minR2 = 20;          % cutoff for vox selection
ROIs = standardROIs;
subjs = {'group'};

% manual set of baseCond + compConds
% [baseCond, compCond], more flexibly defined
comps = [1 2; 1 3; 2 3];
fitSuffix = '_new';

saveFig = 1;

whichStim = 'photo';
whichModel = 'kayCSS';%'cssExpN';%'cssShift';%

hems = {'lh' 'rh'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 11; titleSize = 14;

% now we load in the data from both hemispheres, and threshold across
load(pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems,fitSuffix,task));
ROInum = cellNum(ROIs,info.ROIs);
subjNum = cellNum(subjs,info.subjs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  create supertitle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

titleText = [expt ' (' hemText(hems) ') '];

titleText = [titleText strTogether(subjs) ' (voxels R^2 > ' num2str(minR2) '), ' whichStim ' stim, ' whichModel ' model'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1:
% 1) XY shift magnitude along visual heirarchy
% 2) size change along visual heirarchy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for cc = 1:size(comps,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
niceFig([.1 .1 .8 .55],fontSize);
numPlots = [1 2]; pl = 1;

baseCond = comps(cc,1);
c = comps(cc,2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1) shifts by heirarchy
    subplot(numPlots(1),numPlots(2),pl)
    roiCh = []; roiGroup = [];
    colors = [];
    for r = 1:length(ROIs)
        if length(subjNum) == 1 prfs = subj(subjNum).roi(ROInum(r)); else prfs = roi(ROInum(r)); end
        [~, changes,signs] = plotXYshift2(prfs.fits(baseCond).vox,prfs.fits(c).vox,prfs.fits(1).ppd,prfs.fits(1).res,1);
        roiCh = [roiCh (changes/prfs.fits(1).ppd).*signs];
        roiGroup = [roiGroup repmat({ROIs{r}},1,length(changes))];
        colors = [colors;condColors(r)];
    end
    
    boxplot(roiCh,roiGroup,'colors',colors); hold on; hline(0,'k:');
    set(gca,'box','off','color','none');
    title('XY Shifts (dva)','fontSize',titleSize);
    
    ylabel([prfs.fits(baseCond).cond ' to ' prfs.fits(c).cond],'FontSize',titleSize+8,'FontWeight','bold');
    set(get(gca,'YLabel'),'Visible','on');
    
    pl = pl+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2) size change by heirarchy
    subplot(numPlots(1),numPlots(2),pl)
    roiCh = []; roiGroup = [];
    for r = 1:length(ROIs)
        if length(subjNum) == 1 prfs = subj(subjNum).roi(ROInum(r)); else prfs = roi(ROInum(r)); end
        [changes] = plotSizeChange(prfs.fits(baseCond).vox,prfs.fits(c).vox,prfs.fits(1).ppd,prfs.fits(1).res,1);
        roiCh = [roiCh changes/prfs.fits(1).ppd];
        roiGroup = [roiGroup repmat({ROIs{r}},1,length(changes))];
    end
    
    boxplot(roiCh,roiGroup,'colors',colors); hold on; hline(0,'k:');
    set(gca,'box','off','color','none');
    title('Size Changes (dva)','fontSize',titleSize);
    
    ylabel([prfs.fits(baseCond).cond ' to ' prfs.fits(c).cond],'FontSize',titleSize+8,'FontWeight','bold');
    set(get(gca,'YLabel'),'Visible','on');
    
    pl = pl+1;


superTitle([titleText],titleSize,.97);

if saveFig == 1
    txt = ['acrossROIs_' roi(r).fits(baseCond).cond 'v' roi(r).fits(c).cond];
    if length(hems) == 1
        txt = [txt hems{1} '_']; end
        txt = [whichModel '_' txt];
    niceSave([dirOf(pwd) 'figures/' expt '/acrossROIs/' ],txt,[],info.subjs); % just save pngs, since these can be generated pretty quickly
end
end
