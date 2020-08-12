% loads & plots distributions of XY changes/anything else per each subject

clear all; close all;

subjs = prfSubjs;%{'MG' 'JG' 'TH' 'EM' 'DF' 'SP'};%
expt = 'fixPRF';

saveFig =0;

minR2 = 'r2-20';          % cutoff for vox selection
ROIs= ['V1' standardROIs('face')];%standardROIs;%

whichStim = 'outline';%'photo';%'internal';%'eyes';%
whichModel = 'kayCSS';%'inflipCSSn';%'flipCSSn';%'cssExpN';%'cssShift';%
whichPlot = 'shaded'; %'shaded' or 'bars'

hems = {'rh' 'lh'};
fitSuffix = '';%'_orig';%
yl = [0 .35]; binWidth = .5; % histogram bins

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 11; titleSize = 14;


% now we load in the data from both hemispheres, and threshold across
load(pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems,fitSuffix));
niceFig([.1 .1 .9 .5]);

for r = 1:length(ROIs)
    ROInum = cellNum(ROIs{r},info.ROIs);
    for c = 1:2
        edges = [.25:binWidth:5.25];
        counts = [];
        for s = 1:length(subjs)
            counts(s,:) = histcounts([subj(s).roi(ROInum).fits(c).vox.eccen],edges,'Normalization','probability');
        end
        allCounts{r,c} = counts;
        switch whichPlot
            case 'shaded'
                subplot(1,length(ROIs),r)
                hold on; plot([edges(1:end-1)+binWidth/2],mean(counts),'Color',roiColors(ROIs{r})*(c*.5));
                hold on; e(c) = errorbar([edges(1:end-1)+binWidth/2],mean(counts),se(counts),'v',roiColors(ROIs{r})*(c*.5)); set(e(c),'FaceAlpha',.5);
            case 'bars'
                %subplot(2,length(ROIs),length(ROIs)*(c-1)+r);
                subplot(1,length(ROIs),r)
                if c == 1; mult = .7; else mult = .3; end
                niceBars_colValue([edges(1:end-1)+(mult*binWidth)],mean(counts),se(counts),yl(1),yl(2),colormap('parula')); hold on;
                
        end
    end
        xticks([edges(2:end)-binWidth/2]); xlabel('bin centers');%xticklabels(edges(1:end-1));
        set(gca,'TickDir','out'); ylim(yl); axis square; set(gca, 'box','off');
        title([ROIs{r}]); xlabel('Eccen'); ylabel('Proportion of pRF centers'); %legend(e,{roi(1).fits(1).cond})
    
end


titleText = [expt ' (' hemText(hems) '), voxels R^2 > ' num2str(minR2) '), ' whichStim ' stim, ' whichModel ' model'];


superTitle(titleText,titleSize,.025);

if saveFig == 1
    txt = ['acrossROIs_' whichPlot '_' hemText(hems)];
    niceSave([dirOf(pwd) 'figures/' expt '/centersByEccen/'],txt,[],[],{'svg' 'png'}); % just save pngs, since these can be generated pretty quickly
end

if onLaptop playSound; end