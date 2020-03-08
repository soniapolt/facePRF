% loads & plots model properties (like size x eccen) along the visual heirarchy for cssFit properties
% 7/27/18: changes from lsline to fitline2derror

clear all; close all;


expt = 'fixPRF';

minR2 = 20;          % cutoff for vox selection
ROIs= standardROIs(7);%[standardROIs('EVC') standardROIs('face')];%('face');% 'V3' 'hV4' 'IOG_faces' 'pFus_faces','mFus_faces'
saveFig = 0;
sampleVox = 300;

whichStim = 'outline';
whichModel = 'kayCSS';%'cssExpN';%
fitSuffix = '';

hems = {'rh' 'lh'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 11; titleSize = 14;

load(pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems,fitSuffix));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1: vox scatter + fit lines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f(2) = niceFig([.1 .1 1 .4],fontSize,1);
for r = 1:length(ROIs)
    f(1) = niceFig([.1 .1 .8 .8],fontSize,1);
    
    
    ROInum = cellNum(ROIs{r},info.ROIs);
    subjFit{2} = [];
    for s = 1:length(subj)
        fits = subj(s).roi(ROInum).fits;
        
        for c = 1:length(fits)
                 if c == 1 mult = .25; else mult = 1; end

        figure(f(1)); subplot(2,length(info.subjs)/2,s);
        
        x = [fits(c).vox.eccen]';
        y = [fits(c).vox.size]';
        X = [x ones(length(fits(c).vox),1)];
        
        if ~isempty(X)
        [h,R2] = fitl1line(X,y);
        subjFit{c} = [subjFit{c}; h];
        h1 = plot(x,X*h','Color',condColors(s,1)*mult); hold on;
        hold on; scatter(x,y,5,condColors(s,1)*mult,'filled'); hold on; title([ROIs{r} ' subj ' info.subjs{s}]);
        % figure(f(2)); subplot(1,length(ROIs),r); sc{s} = scatter(x,y,5,roiColors(ROIs{r})*mult,'filled'); hold on; 
        end
        end
        figure(f(1));
    xl = xlim; xlim([0.25 6]); yl = ylim; if containsTxt(ROIs{r},'faces') ylim([0 18]); else ylim([0 12]);end
    set(gca,'TickDir','out');
    xlabel('Eccen (dva)','FontSize',fontSize); ylabel(['Size (2*SD/sqrt(N)) (dva)'],'FontSize',fontSize);
    axis square;
    
    end
    
    
    
    figure(f(2));subplot(1,length(ROIs),r);
    fits = roi(ROInum).fits;
        
    
    if sampleVox < 1 sv = round(length(fits(1).vox)*sampleVox);
        elseif sampleVox > length(fits(1).vox) sv = length(fits(1).vox);
        else sv = sampleVox; end
        sv = datasample([1:length(fits(1).vox)],sv,'Replace',false);
    
    for c = 1:length(fits)
        if c == 1 mult = .25; else mult = 1; end
            s(c) = scatter([fits(c).vox(sv).eccen],[fits(c).vox(sv).size],2,roiColors(ROIs{r})*mult,'filled'); hold on;

      meanH = mean(vertcat(subjFit{c}));
      seH = se(vertcat(subjFit{c}));
      x = [0:.1:6]'; X = [x ones(length(x),1)];
      h1 = plot(x,X*meanH','Color',roiColors(ROIs{r})*mult); hold on;
       h2 = plot(x,X*(meanH+seH)'); set(h2,'Color',roiColors(ROIs{r})*mult,'LineStyle',':'); hold on;
       h3 = plot(x,X*(meanH-seH)'); set(h3,'Color',roiColors(ROIs{r})*mult,'LineStyle',':'); hold on;

    end
    
    xl = xlim; xlim([0.25 6]); yl = ylim; if containsTxt(ROIs{r},'faces') ylim([0 18]); else ylim([0 12]);end
    set(gca,'TickDir','out');
    xlabel('Eccen (dva)','FontSize',fontSize); ylabel(['Size (2*SD/sqrt(N)) (dva)'],'FontSize',fontSize);
    axis square; %legend(sc{:},prfSubjs);
    title(ROIs{r});

if saveFig == 1
    txt = [ROIs{r} '_' whichModel '_' whichStim '_r' num2str(minR2) '_' hemText(hems)]; 
    figure(f(1)); superTitle(txt,12,.97);
    niceSave([dirOf(pwd) 'figures/' expt '/sizeEccen/randomE/'],[ROIs{r} '_' txt],[],[],{'png' 'svg'}); % just save pngs, since these can be generated pretty quickly
end
end
    if saveFig == 1 figure(f(2));
    niceSave([dirOf(pwd) 'figures/' expt '/sizeEccen/randomE/'],['acrossSubjs_' txt],[],[],{'png' 'svg'}); % just save pngs, since these can be generated pretty quickly
    end
playSound;