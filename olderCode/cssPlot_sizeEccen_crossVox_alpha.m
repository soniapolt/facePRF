% loads & plots model properties (like size x eccen) along the visual heirarchy for cssFit properties
% 7/27/18: changes from lsline to fitline2derror

clear all; close all;


expt = 'fixPRF';

minR2 = 20;          % cutoff for vox selection
ROIs= standardROIs(7);%['V1' standardROIs('face')];%('face');% 'V3' 'hV4' 'IOG_faces' 'pFus_faces','mFus_faces'
sampleVox = 300; 
saveFig = 0;

whichStim = 'outline';
whichModel = 'kayCSS';%'cssExpN';%
fitSuffix = '';

hems = {'rh' 'lh'};
yl = 12;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 11; titleSize = 14;

load(pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems,fitSuffix));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f(1) = niceFig([.1 .1 .4 .4],fontSize,1);
f(2) = niceFig([.1 .1 .7 .4],fontSize,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1: vox scatter + fit lines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for r = 1:length(ROIs)
    figure(f(1));
    subplot(2,ceil(length(ROIs)/2),r);
    
    ROInum = cellNum(ROIs{r},info.ROIs);
    fits = roi(ROInum).fits;
    
    
    if sampleVox < 1 sv = round(length(fits(1).vox)*sampleVox);
        elseif sampleVox > length(fits(1).vox) sv = length(fits(1).vox);
        else sv = sampleVox; end
        sv = datasample([1:length(fits(1).vox)],sv,'Replace',false);

    for c = 1:length(fits)
        if c == 1 mult = .15; else mult = 1; end
    alpha = mat2gray([fits(c).vox(sv).r2]); %ones(1,length(fits(c).vox));% 
    hold on;
    s{c} = scatterAlpha([fits(c).vox(sv).eccen],[fits(c).vox(sv).size]/2,alpha,roiColors(ROIs{r})*mult,2); hold on;
    
    fitRange = [0:.5:8];
    
    
                x = [fits(c).vox.eccen]';
                y = [fits(c).vox.size]';
                X = [x ones(length(fits(c).vox),1)];
                
                if ~isempty(X)
                    [h,R2] = fitl1line(X,y);
                    h1 = plot(x,X*h','Color',roiColors(ROIs{r})*mult); hold on;
                    s{c} = scatterAlpha([fits(c).vox(sv).eccen],[fits(c).vox(sv).size]/2,alpha,roiColors(ROIs{r})*mult,2); hold on;
                    set(h1,'LineWidth',1); %set(extLine,'LineStyle',':');
                    %alpha(.2);
                end
    
    %if containsTxt(ROIs{r},'faces') fitRange = [.25:.25:4]; else fitRange = [.25:.25:6];end
% % %     [N,edges] = histcounts([fits(c).vox.eccen],fitRange);
% % %     [errorObj,lineObj] = scatterline([fits(c).vox.eccen],[fits(c).vox.size],fitRange(N>9),NaN,100,roiColors(ROIs{r})*mult,2,1); hold on;
% % %     extX = [0:.1:6];
% % %     extLine = plot(extX, extendLine(lineObj,extX),'Color',roiColors(ROIs{r})*mult,'LineWidth',.5);
% % %     set(lineObj,'LineWidth',.5); set(extLine,'LineStyle',':');
%     alpha(.2);
    end
    xlim([0 6]); ylim([0 6]);
    %xl = xlim; xlim([0 6]); y = ylim; if containsTxt(ROIs{r},'faces') ylim([0 6]); else ylim([0 6]);end
    set(gca,'TickDir','out');
    xlabel('Eccen (dva)','FontSize',fontSize); ylabel(['Size (SD/sqrt(N)) (dva)'],'FontSize',fontSize); 
    axis square; 
    title(ROIs{r});
    
end
%line; l=fliplr(l); for r = 1:length(ROIs) set(l(r),'Color', condColors(r,1),'LineWidth',2); end

%xlim([0 12]); ylim([0 12]);


%superTitle(titleText,titleSize);

if saveFig == 1
    txt = ['sizeEccen_' whichModel '_' whichStim '_r' num2str(minR2) '_' hemText(hems)];
    figure(f(1));
    niceSave([dirOf(pwd) 'figures/' expt '/sizeEccen/'],['scatter_' txt '_ylim' num2str(yl)],[],[],{'png' 'svg'}); % just save pngs, since these can be generated pretty quickly
    figure(f(2));
    niceSave([dirOf(pwd) 'figures/' expt '/sizeEccen/'],['lines_' txt '_ylim' num2str(yl)],[],[],{'png' 'svg'}); % just save pngs, since these can be generated pretty quickly
end
 playSound;