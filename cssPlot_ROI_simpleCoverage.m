% loads & plots some basics re: invPRF fits (css model)
% SP 10/25/17
% _indiv code allows us to plot individual ROIs in multiple plots
% 10/8/18 edit: made more flexible across expriment with different numbers
% of fit conditions

clear all; close all;

subjs = {'SP' 'DF' 'EM' 'TH' 'MG' 'JG'};%;
expt = 'fixPRF';

minR2 = 20;          % cutoff for vox selection
ROIs= standardROIs; %
sampleVox = 0; % how many randomly selected voxels are we plotting?

saveFig = 1;

whichStim = 'photo';%'internal';%
whichModel = 'kayCSS';%''cssExpN';%cssShift';%
fitSuffix = '';

hems = {'lh' 'rh'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 12; titleSize = 14;

% now we load in the data from both hemispheres, and threshold across
load(pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems,fitSuffix));
ROInum = cellNum(ROIs,info.ROIs);
subjNum = cellNum(subjs,info.subjs);

for r = 1:length(ROIs)
    if length(subjNum) == 1 bFits = subj(subjNum).roi(ROInum(r)).fits;
    else bFits =roi(ROInum(r)).fits; end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  create supertitle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if exist('task','var');
        titleText = [task ' task, ']; else titleText = [];end
    titleText = [whichModel ' ' titleText hemText(hems) ' ' ROIs{r} ', Subjs: ' strTogether(subjs) ...
        ' (' num2str(length(bFits(1).vox)) ' voxels R^2 > ' num2str(minR2) ', sampleVox = ' num2str(sampleVox) '), ' whichStim];
    
    if onLaptop niceFig([.1 .1 .8 .8],fontSize); else
        niceFig([0 .2 .6 .25*(size(comps,1)+1)],fontSize); end
    numPlots = [1 length(bFits)]; pl = 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1) coverage, 
    for c = fliplr([1:length(bFits)])
    subplot(numPlots(1),numPlots(2),pl)
    plotCoverage(bFits(c).vox,condColors(c),'',roi(1).fits(1).ppd,roi(1).fits(1).res,sampleVox,1);
    t = title(roi(1).fits(c).cond);
    set(t,'visible','on');
    pl = pl+1;
    end
    
    superTitle(titleText,titleSize,.05);
    if saveFig
        txt = [ROIs{r} '_' whichModel  '_' whichStim '_sampleVox' num2str(sampleVox)];
        niceSave([dirOf(pwd) 'figures/' expt '/simpleCoverage/'],txt,[],subjs); % just save pngs, since these can be generated pretty quickly
    end
    
    
end % ROIs

if onLaptop playSound; end