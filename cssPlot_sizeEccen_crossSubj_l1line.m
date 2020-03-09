% loads & plots model properties (like size x eccen) along the visual heirarchy for cssFit properties
% 7/27/18: changes from lsline to fitline2derror

clear all; close all;


expt = 'fixPRF';

minR2 = 20;          % cutoff for vox selection
ROIs= standardROIs('face');%[standardROIs('EVC') standardROIs('face')];%('face');% 'V3' 'hV4' 'IOG_faces' 'pFus_faces','mFus_faces'
saveFig = 1;
sampleVox = 300;
fitIndivs = 1;

whichStim = 'outline';
whichModel = 'kayCSS';%'cssExpN';%
fitSuffix = '';

hems = {'rh' 'lh'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 11; titleSize = 14;
pf = pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems,fitSuffix);
load(pf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1: vox scatter + fit lines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f(2) = niceFig([.1 .1 1 .4],fontSize,1);
for r = 1:length(ROIs)
    outFile = ['sizeEccen/' hemText(hems) '_' ROIs{r} '_' whichStim '_' whichModel '_r2-' num2str(minR2) '.mat'];
    if ~exist(outFile) || fitIndivs == 1
        f(1) = niceFig([.1 .1 .8 .8],fontSize,1);
        
        ROInum = cellNum(ROIs{r},info.ROIs);
        sFit = struct;
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
                    h1 = plot(x,X*h','Color',condColors(s,1)*mult); hold on;
                    hold on; scatter(x,y,5,condColors(s,1)*mult,'filled'); hold on; title([ROIs{r} ' subj ' info.subjs{s}]);
                    set(h1,'LineWidth',1); %set(extLine,'LineStyle',':');
                    alpha(.2);
                    
                    % save this fit
                    sFit(c).subj{s} = sort(horzcat(x,X*h'));
                    
                end
            end
            figure(f(1));
            xl = xlim; xlim([0.25 6]); yl = ylim; if containsTxt(ROIs{r},'faces') ylim([0 18]); else ylim([0 12]);end
            set(gca,'TickDir','out');
            xlabel('Eccen (dva)','FontSize',fontSize); ylabel(['Size (2*SD/sqrt(N)) (dva)'],'FontSize',fontSize);
            axis square;
            
        end
        
        save(outFile,'sFit','pf'); fprintf('Saving %s...\n',outFile);
        txt = ['fitl1lineversion_' ROIs{r} '_' whichModel '_' whichStim '_r' num2str(minR2) '_' hemText(hems)];
        if saveFig == 1
            figure(f(1)); superTitle(txt,12,.97);
            niceSave([dirOf(pwd) 'figures/' expt '/sizeEccen/randomE/'],txt,[],[],{'png' 'svg' 'fig'}); % just save pngs, since these can be generated pretty quickly
        end
    else
        load(outFile); fprintf('Loading %s...\n',outFile);
        txt = ['fitl1lineversion_' ROIs{r} '_' whichModel '_' whichStim '_r' num2str(minR2) '_' hemText(hems)];
        try load([dirOf(pwd) 'figures/' expt '/sizeEccen/randomE/' txt '.fig'])
        catch warning('Unable to load a .fig for these subject fits!');
        end
    end
    
    figure(f(2)); subplot(1,length(ROIs),r);
    for c = 1:2
        if c == 1 mult = .25; else mult = 1; end
        a = vertcat(sFit(c).subj{:});
        x = []; y = [];
        for s = 1:12
            b = sFit(1).subj{s};
            x(end+1) = b(1,1); y(end+1) = b(1,2);
            x(end+1) = b(end,1); y(end+1) = b(end,2);
        end
        x = a(:,1); y = a(:,2);
        
        
        fitRange = [0:.5:8];
        [N,edges] = histcounts(y,fitRange);
        scatter(x,y,2,roiColors(ROIs{r})*mult,'filled'); hold on;
        [errorObj,lineObj,mn,se] = scatterline(x,y,fitRange(N>1),NaN,1000,roiColors(ROIs{r})*mult,2,1); hold on;
        %delete(errorObj);
        extX = [0:.1:8];
        extLine = plot(extX, extendLine(lineObj,extX),'Color',roiColors(ROIs{r})*mult,'LineWidth',.5);
        set(lineObj,'LineWidth',.5); set(extLine,'LineStyle',':');
        alpha(.2);
        title(ROIs{r});
    end
end
playSound;