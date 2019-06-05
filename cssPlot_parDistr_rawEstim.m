% loads & plots distributions of XY changes/anything else per each subject

clear all; close all;

subjs = {'MG' 'JG' 'TH' 'EM' 'DF' 'SP'};%
task = '';
expt = 'fixPRF';
noCenters = 0;

saveFig = 0;
convertDVA = 1;
normVoxCount = 1;

minR2 = 20;          % cutoff for vox selection
ROIs= {'mFus_faces'};%'V1' 'V2' 'V3' 'hV4' hV4' 'IOG_faces' 'pFus_faces' 

whichStim = 'photo';%'eyes';%'internal';%
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
% figure 1: distribution of parameters for this ROI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p = 1:length(roi(1).fits(1).parNames)
    
titleText = [whichModel ' ' roi(1).fits(1).parNames{p} ', Subj: '];
titleText = [titleText strTogether(subjs) ' (voxels R^2 > ' num2str(minR2) '), ' whichStim ' stim, ' whichModel ' model'];
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if saveFig && onLaptop figSize = [0 0 1 1]; else figSize = [.2 .1 8 .8]; end
    niceFig(figSize,fontSize,1);
    numPlots = [2 ceil(length(ROIs)/2)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for r = 1:length(ROIs)
    subplot(numPlots(1),numPlots(2),r)
    
    
    for c = 1:length(roi(1).fits)
        pars = vertcat(roi(ROInum(r)).fits(c).vox.params);   
        % get param values for this condition
        
        cPars{c} = pars(:,p)';
        if convertDVA && ~containsTxt(roi(1).fits(1).parNames{p},'exp') && ~containsTxt(roi(1).fits(1).parNames{p},'gain')
        % rescale some parameters so that they are in DVA units and
        % centered around zero (center of screen)
        if ~containsTxt(roi(1).fits(1).parNames{p},'sd') % don't re-center the SD
            cPars{c} = roi(1).fits(1).res-cPars{c}-roi(1).fits(1).res/2;
        end
        cPars{c} = cPars{c}./roi(1).fits(1).ppd;
        end
    end   
    if containsTxt(roi(1).fits(1).parNames{p},'gain')
        niceBoxplot([cPars{1};cPars{2}]',{roi(ROInum(r)).fits.cond},1,[condColors(4);condColors(2)],[0 10]);
        %boxplot([cPars{1};cPars{2}]','labels',{roi(ROInum(r)).fits.cond},'plotstyle','traditional','outliersize',.5,'jitter',.1);
        %ylim([0 10]);
    else
        
    plotDistr(cPars,1,{roi(ROInum(r)).fits.cond},nBins,whichM,0);
    end
    %title(ROIs{r},'fontSize',titleSize,'interpreter','none','FontWeight','bold'); 
    %xlabel(roi(1).fits(1).parNames{p},'fontSize',titleSize);
    xlabel(ROIs{r},'fontSize',titleSize,'interpreter','none','FontWeight','bold');

    end
    superTitle(titleText,titleSize,.025);

if saveFig == 1
    if length(subjs) == 1 txt = ['subj' subjs{1}]; else txt = ['groupN' num2str(length(subjs))]; end
    switch expt
            case 'invPRF3'
                txt = [roi(1).fits(1).parNames{p} '_' task 'Task_' txt ];
            case 'fixPRF'
                txt = [roi(1).fits(1).parNames{p} '_' txt ];
            case 'compPRF'
                txt = [roi(1).fits(1).parNames{p} '_' txt ];
    end
    if length(hems) == 1
        txt = [txt '_' hems{1}]; end
        txt = ['distr_' whichModel '_' whichStim '_' txt];
    niceSave([dirOf(pwd) 'figures/' expt '/params/'],txt); % just save pngs, since these can be generated pretty quickly
end
end
if onLaptop playSound; end