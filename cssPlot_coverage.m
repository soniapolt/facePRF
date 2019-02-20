% loads & plots changes along the visual heirarchy for cssFit properties

clear all; close all;

subjs = {'JG' 'MG' 'DF' 'TH' 'EM' 'SP'};%
tasks = {'fix'};
expt = 'fixPRF';
fitSuffix = '';

saveFig = 1;

minR2 = 20;          % cutoff for vox selection
ROIs= standardROIs;%{'V1' 'hV4' 'IOG_faces' 'mFus_faces' 'pFus_faces'};
conds = [1 2]; % 1 = inverted, 2 = misaligned, 3 = normal, or 1 = inverted 2 = misaligned
covLims = [0 1]; % min and max for plotting coverage

whichStim = 'internal';%'photo';%
whichModel = 'kayCSS';%'cssShift';%
hems = {'rh' 'lh'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load in data, sort plotVox           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 11; titleSize = 14;

for t = 1:length(tasks)
% now we load in the data from both hemispheres, and threshold across
load(pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems,fitSuffix));
task(t).roi = roi;
end % tasks

ROInum = cellNum(ROIs,info.ROIs);
subjNum = cellNum(subjs,info.subjs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  create supertitle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

titleText = [expt ' ' strTogether(subjs) ' (' hemText(hems) ' voxels R^2 > ' num2str(minR2) '), ' whichStim ' stim, ' whichModel ' model'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if onLaptop figSize = [0 0 1 .33*length(conds)];
else figSize = [.2 .2 .5 .15*length(conds)];
end
niceFig(figSize,fontSize);
numPlots = [length(conds) length(ROIs)*length(tasks)];  pl = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1:
% rows are conditions
% colums are ROI*task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for cc = 1:length(conds)
    c  = conds(cc);
    for r = 1:length(ROIs)
        
        for t = 1:length(tasks)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % overall coverage
            subplot(numPlots(1),numPlots(2),pl); pl = pl+1;
            
            sumCoverage(task(t).roi(ROInum(r)).fits(c).vox,roi(1).fits(1).ppd,roi(1).fits(1).res,covLims,1);
            
            set(get(gca,'Title'),'Visible','on');
            title({ROIs{r}; [tasks{t} ' task']},'fontSize',titleSize,'Interpreter','none','Color',condColors(r+3,1));
            % t = title(['{\color[rgb]{' num2str(condColors(r+3,1)) '}' ROIs{r} '} ' tasks{t} ' task']);
            
            if r == 1 && t == 1 % left edge of plots
                set(get(gca,'YLabel'),'Visible','on');
                ylabel(task(t).roi(r).fits(c).cond,'FontSize',titleSize+8,'FontWeight','bold');
            end
        end
    end
end

superTitle([titleText],titleSize,.05);

if saveFig == 1
    txt = [whichModel '_' 'coverage_' roi(1).fits(conds(1)).cond 'v' roi(1).fits(conds(2)).cond];
    if ~containsTxt(whichStim,'photo')
        txt = [ txt '_' whichStim ]; end
    niceSave([dirOf(pwd) 'figures/' expt '/coverage/'],txt,ROIs,subjs); % just save pngs, since these can be generated pretty quickly
end

%playSound;
