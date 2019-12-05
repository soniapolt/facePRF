% plots parameters for each ROI (as apposed to across ROIs)

clear all; close all;

subjs = {'TH' 'DF' 'EM' 'MG' 'JG' 'SP'};%
task = '';
expt = 'fixPRF';

saveFig = 1;

minR2 = 20;          % cutoff for vox selection
ROIs= standardROIs;%{'IOG_faces' 'pFus_faces' 'mFus_faces'};%'V1' 'V2' 'V3' 'hV4' 

whichStim = 'internal';%'photo';
whichModel = 'kayCSS';%'cssShift';%
whichM = 3; % 1 = mean, 2 = mode/peak, 3 = median
hems = {'rh' 'lh'};
fitSuffix = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 11; titleSize = 14;
nBins = 20; % histogram bins
load(pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems,fitSuffix,task));
ROInum = cellNum(ROIs,info.ROIs);
subjNum = cellNum(subjs,info.subjs);

if length(subjNum) == 1 roi = subj(subjNum).roi; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1: distribution of parameters for this ROI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for r = 1:length(ROIs)
titleText = [hemText(hems) ' ' ROIs{r} ', Subj: '];
titleText = [titleText strTogether(subjs) ' (voxels R^2 > ' num2str(minR2) '), ' whichStim ' stim, ' whichModel ' model'];
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if saveFig && onLaptop figSize = [0 0 1 1]; else figSize = [.2 .1 .8 .8]; end
    niceFig(figSize,fontSize);
    numPlots = [2 ceil(length(roi(1).fits(1).parNames)/2)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1) shifts across ROIs
    for p = 1:length(roi(1).fits(1).parNames)
    subplot(numPlots(1),numPlots(2),p)
    
    for c = 1:length(roi(1).fits)
        pars = vertcat(roi(ROInum(r)).fits(c).vox.params);   
        % get param values for this condition
        cPars{c} = pars(:,p)';
    end
    
    plotDistr(cPars,1,{roi(ROInum(r)).fits.cond},nBins,whichM);
    
    title(roi(1).fits(1).parNames{p},'fontSize',titleSize,'interpreter','none','FontWeight','bold'); 
    
    end
    superTitle(titleText,titleSize,.97);

if saveFig == 1
    switch expt
            case 'invPRF3'
                txt = ['pars_' task 'Task_' ROIs{r}];
            case 'fixPRF'
                txt = ['pars_' ROIs{r}];
            case 'compPRF'
                txt = ['pars_' ROIs{r}];
    end
    if length(hems) == 1
        txt = [txt '_' hems{1}]; end
    
    if ~containsTxt(whichStim,'photo')
        txt = [whichStim '_' txt]; end
        
        txt = [whichModel '_' txt];
    niceSave([dirOf(pwd) 'figures/' expt '/paramsROI/'],txt,[],subjs); % just save pngs, since these can be generated pretty quickly
end
end