% loads & plots changes along the visual heirarchy for cssFit properties
% this version is for performing stats, not plotting
clear all; close all;

subjs = {'TH' 'JG' 'DF' 'MG' 'EM' 'SP'};
tasks = {'fix'};
expt = 'fixPRF';

saveFig = 1;

minR2 = 20;          % cutoff for vox selection
ROIs= standardROIs;
conds = [2 1]; 
covLims = [0 1]; % min and max for plotting coverage
fitSuffix = '';%'_orig';%

whichStim = 'photo';%'internal';%
whichModel = 'kayCSS';%'cssExpN';%'cssShift';%

hems = {'lh'};

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
titleText = [ whichModel ' ' hemText(hems) '_ ' ROIs{r} ' (R^2 > ' num2str(minR2) '), ' whichStim ' stim, ' whichModel ' model'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if onLaptop figSize = [0 0 1 .33*length(conds)];
else figSize = [.2 .2 .5 .15*length(conds)];
end
niceFig(figSize,fontSize);
numPlots = [length(conds) length(subjs)];  pl = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1:
% rows are conditions
% colums are subj*task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for cc = 1:length(conds)
    c  = conds(cc);
    
        for s = 1:length(subjs)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % overall coverage
            subplot(numPlots(1),numPlots(2),pl); pl = pl+1;
            prfs = subj(subjNum(s)).roi(ROInum(r)); 
            [coverage, areaDeg(s,cc)] = sumCoverage(prfs.fits(c).vox,prfs.fits(1).ppd,prfs.fits(1).res,covLims,1);
            
            set(get(gca,'Title'),'Visible','on');
            title({subjs{s}; [num2str(areaDeg(s,cc)) ' deg^2']},'fontSize',titleSize,'Interpreter','none','Color',condColors(r+3,1));
            % t = title(['{\color[rgb]{' num2str(condColors(r+3,1)) '}' ROIs{r} '} ' tasks{t} ' task']);
            
            if s == 1  % left edge of plots
                set(get(gca,'YLabel'),'Visible','on');
                ylabel(prfs.fits(c).cond,'FontSize',titleSize+8,'FontWeight','bold');
            end
        end
end
superTitle([titleText],titleSize,.05);

% stats on coverage differences
[H,P,CI,STATS] = ttest(areaDeg(:,1),areaDeg(:,2))

superTitle(['t(' num2str(length(subjs)-1) ')=' num2str(STATS.tstat) ', p=' num2str(P)],titleSize,.01);

if saveFig == 1
    txt = [whichModel]; 
    if ~containsTxt(whichStim,'photo')
        txt = [ txt '_' whichStim ]; end
    if length(hems) == 1
        txt = [txt '_' hems{1}];
    end
    txt = [txt '_' prfs.fits(conds(1)).cond 'v' prfs.fits(conds(2)).cond '_'];
    niceSave([dirOf(pwd) 'figures/' expt '/coverage2/'],txt,ROIs(r),subjs); % just save pngs, since these can be generated pretty quickly
end
end

if onLaptop playSound;end
