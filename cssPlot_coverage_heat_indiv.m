% loads & plots changes along the visual heirarchy for cssFit properties
% this version is for performing stats, not plotting
clear all; close all;

subjs = prfSubjs;%{ 'JP' 'MH' 'JW'};%'TH' 'JG' 'DF' 'MG' 'EM' 'SP'
expt = 'fixPRF';

saveFig = 1;

minR2 = 20;          % cutoff for vox selection
ROIs= standardROIs('face+');%8);%
conds = [2 1];
covLims = [0 1]; % min and max for plotting coverage
fitSuffix = '';%'_orig';%

whichStim = 'outline';%'photo';%'internal';%
whichModel = 'kayCSS';%'cssExpN';%'cssShift';%

hems = {'lh' 'rh'};
figDir = [dirOf(pwd) 'figures/' expt '/subjCov/'];

plotIndivs = 1;
plotSummary = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load in data, sort plotVox           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fontSize = 11; titleSize = 14;
% now we load in the data from both hemispheres, and threshold across
load(pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems,fitSuffix));
ROInum = cellNum(ROIs,info.ROIs);
subjNum = cellNum(subjs,info.subjs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  create supertitle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for r = 1:length(ROIs)
    titleText = [whichModel ' ' hemText(hems) '_ ' ROIs{r} ' (R^2 > ' num2str(minR2) '), ' whichStim ' stim, ' whichModel ' model'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if onLaptop figSize = [0 0 1 .33*length(conds)];
    else figSize = [.2 .2 .5 .5];
    end
    % across-subj figure
    if plotSummary
        mf = niceFig([.1 .1 .9 .5],fontSize); numPlots = [length(conds) length(subjs)];  pl = 1; end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % figure 1:
    % rows are conditions
    % colums are subj*task
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for s = 1:length(subjs)
        if plotIndivs of = niceFig(figSize,fontSize);end
        
        for cc = 1:length(conds)
            c  = conds(cc);
            prfs = subj(subjNum(s)).roi(ROInum(r));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % overall coverage
            if plotIndivs
                figure(of); subplot(1,2,c)
                [coverage, areaDeg(s,cc)] = sumCoverage(prfs.fits(c).vox,prfs.fits(1).ppd,prfs.fits(1).res,covLims,1);
                set(get(gca,'Title'),'Visible','on');
                title({[subjs{s} ' ' num2str(length(prfs.fits(c).vox)) ' voxels'];...
                    prfs.fits(c).cond; [num2str(areaDeg(s,cc)) ' deg^2']},'fontSize',titleSize,'Interpreter','none','Color',condColors(r+3,1));
                superTitle([titleText],titleSize,.05);
            end
            
            if plotSummary
                figure(mf);sp = subplot(numPlots(1),numPlots(2),numPlots(2)*(cc-1)+s); 
                [coverage, areaDeg(s,cc)] = sumCoverage(prfs.fits(c).vox,prfs.fits(1).ppd,prfs.fits(1).res,covLims,1);
                
                set(get(gca,'Title'),'Visible','on');
                title({subjs{s};prfs.fits(c).cond; [num2str(areaDeg(s,cc)) ' deg^2']},'fontSize',titleSize,'Interpreter','none','Color',condColors(r+3,1));
                
                if s == 1  % left edge of plots
                    set(get(gca,'YLabel'),'Visible','on');
                    ylabel(prfs.fits(c).cond,'FontSize',titleSize+8,'FontWeight','bold');
                end
                superTitle([titleText],titleSize,.05);
            end
            
        end
        
        % stats on coverage differences
        [H,P,CI,STATS] = ttest(areaDeg(:,1),areaDeg(:,2))
        fprintf(['t(' num2str(length(subjs)-1) ')=' num2str(STATS.tstat) ', p=' num2str(P) '\n'],titleSize,.01);
        
        % base title for figs
            btxt = [whichModel];
            if ~containsTxt(whichStim,'photo')
                btxt = [ btxt '_' whichStim ]; end
            if length(hems) == 1
                btxt = [btxt '_' hems{1}];
            end
            
            if plotIndivs && saveFig
                figure(of)
                txt = [hemText(hems) '_' ROIs{r} '_' subjs{s} '_' btxt '_' prfs.fits(conds(1)).cond 'v' prfs.fits(conds(2)).cond  ];
                niceSave(figDir,txt); % just save pngs, since these can be generated pretty quickly
            end
        end

    
    if plotSummary && saveFig
        figure(mf);
        superTitle(['t(' num2str(length(subjs)-1) ')=' num2str(STATS.tstat) ', p=' num2str(P)],titleSize,.01);
        
        txt = ['group_' ROIs{r} '_' num2str(length(subjs)) '_' btxt '_' prfs.fits(conds(1)).cond 'v' prfs.fits(conds(2)).cond];
        niceSave(figDir,txt); % just save pngs, since these can be generated pretty quickly
    end
    close all;save([figDir 'N' num2str(length(subjs)) '_areaDeg_' ROIs{r} '.mat'],'areaDeg','H','P','STATS');
end   
    

if onLaptop playSound;end
