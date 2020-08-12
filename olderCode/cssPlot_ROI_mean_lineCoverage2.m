% loads & plots some basics re: invPRF fits (css model)
% SP 10/25/17
% _indiv code allows us to plot individual ROIs in multiple plots
% 10/8/18 edit: made more flexible across expriment with different numbers
% of fit conditions

clear all; close all;

subjs = {'SP' 'DF' 'EM' 'TH' 'MG' 'JG'};%;
expt = 'fixPRF';

minR2 = 20;          % cutoff for vox selection
ROIs= standardROIs('face'); %

saveFig = 1;

whichStim = 'photo';%'internal';%
whichModel = 'kayCSS';%''cssExpN';%cssShift';%
whichM = 'median';

fitSuffix = '';

hems = {'rh' 'lh'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 12; titleSize = 14;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  create supertitle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('task','var');
    titleText = [task ' task, ']; else titleText = [];end
titleText = [hemText(hems) ' ' whichModel ' ' titleText ', Subjs: ' strTogether(subjs) ...
    ' (R^2 > ' num2str(minR2) '), ' whichStim];

if onLaptop niceFig([.1 .1 .8 .8],fontSize,1); else
    niceFig([0 .2 .6 .25*(size(comps,1)+1)],fontSize,1); end
numPlots = [1 2]; pl = 1;
% now we load in the data from both hemispheres, and threshold across
load(pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems,fitSuffix));

ROInum = cellNum(ROIs,info.ROIs);
subjNum = cellNum(subjs,info.subjs);

% 1) coverage,
    for c = fliplr([1:length([roi(1).fits])])
        %subplot(numPlots(1),numPlots(2),length(ROIs)*(c-1)+r)
        subplot(numPlots(1),numPlots(2),c)
        plotMeanCoverage2(subj,whichM,ROInum,c,3.2,ROIs);
        t = title([roi(1).fits(c).cond]);
        set(t,'visible','on');
    end
    
    %subplotresize(numPlots(1), numPlots(2),.7,.7);
    
    superTitle(titleText,titleSize,.05);
    if saveFig
        
        txt = ['mean_' hemText(hems) '_' strTogether(ROIs) '_' whichModel  '_' whichStim];
        if ~isequal(minR2,20)
            txt = [txt '_r2-' num2str(minR2)];end
        niceSave([dirOf(pwd) 'figures/' expt '/lineCoverage/'],txt,[],subjs); % just save pngs, since these can be generated pretty quickly
    end
    
if onLaptop playSound; end