% loads & plots model properties (like size x eccen) along the visual heirarchy for cssFit properties
% 7/27/18: change from lsline to fitline2derror
% 3/10/20: change from fitline2derror/scatterline to bootstrapping over
% fitl1line
% 5/6/20: size as 1 sigma/sqrt(n) change
%%% THIS IS THE FIG 2 CODE

clear all; close all;


expt = 'fixPRF';

minR2 = 50;          % cutoff for vox selection
ROI= standardROIs(7);%('face');% 'V3' 'hV4' 'IOG_faces' 'pFus_faces','mFus_faces'
scatterVox = 300;
doAlpha = 0;
saveFig = 1;
bootMethod = '1l';%'scatterline';%'1l'

whichStim = 'outline';
whichModel = 'kayCSS';%'cssExpN';%
fitSuffix = '';

hems = {'rh' 'lh'};
yl = 6;

boot.numIter = 100;
boot.sampleVox = .8;
boot.binThresh = 5; % need n+1 voxels to bootstrap errorbar
boot.fitRange = [0:.5:10];%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 11; titleSize = 14;

load(pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems,fitSuffix));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%f(1) = niceFig([.1 .1 .8 .4],fontSize,1);
f(1) = niceFig([.1 .1 .25 .4],fontSize,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1: vox scatter + fit lines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    ROInum = cellNum(ROI,info.ROIs);
    fits = roi(ROInum).fits;
    
    if scatterVox < 1 sv = round(length(fits(1).vox)*scatterVox);
    elseif scatterVox > length(fits(1).vox) sv = length(fits(1).vox);
    else sv = scatterVox; end
    sv = datasample([1:length(fits(1).vox)],sv,'Replace',false);
    
    for c = 1:length(fits)
        figure(f(1)); %subplot(2,round(length(ROIs)/2),r);
        if c == 1 mult = .15; else mult = 1; end
        if doAlpha
            alpha = mat2gray([fits(c).vox(sv).r2]);
        else alpha = ones(1,length(sv));end
        hold on;
        s{c} = scatterAlpha([fits(c).vox(sv).eccen],[fits(c).vox(sv).size]/2,alpha,roiColors(ROI)*mult,4); hold on;
        
        
        %%% bootstrapping stage
        switch bootMethod
            case '1l'
                
                h = NaN(boot.numIter,2); R2 = NaN(boot.numIter,1);
                
                % points at which we'll compute the bootstrapped std
                [N,edges] = histcounts([fits(c).vox.eccen],boot.fitRange);
                fitRange = [boot.fitRange(N>boot.binThresh)' ...
                    ones(length(boot.fitRange(N>boot.binThresh)),1)]; fitOut = NaN(boot.numIter,length(fitRange));
                
                %tic
                parfor b = 1:boot.numIter
                    v = datasample([1:length(fits(1).vox)],round(length(fits(1).vox)*boot.sampleVox),'Replace',true);
                    x = [fits(c).vox(v).eccen]';
                    y = [fits(c).vox(v).size]'/2;
                    X = [x ones(length(v),1)];
                    [h(b,:),R2(b)] = fitl1line(X,y);
                    
                    fitOut(b,:) = fitRange*h(b,:)';
                end
                %toc
                
                % full line fit
                x = [fits(c).vox.eccen]';
                y = [fits(c).vox.size]'/2;
                X = [x ones(length(fits(c).vox),1)];
                
                [hFull,R2Full] = fitl1line(X,y);
                h1 = plot(boot.fitRange,[boot.fitRange' ones(length(boot.fitRange),1)]*hFull','Color',roiColors(ROI)*c*.5); hold on;
                [ci,med] = CI(fitOut);
                
                % construct errorbar
                plotErr = NaN(2,length(boot.fitRange));
                plotErr(:,find(N>boot.binThresh)) = ci;
                
                hold on; e = errorbar3(boot.fitRange,[[boot.fitRange' ones(length(boot.fitRange),1)]*hFull']',plotErr,'v',roiColors(ROI)*c*.5); set(e,'FaceAlpha',.2);
                
            case 'scatterline'
                x = [fits(c).vox.eccen]';
                y = [fits(c).vox.size]'/2;
                %X = [x ones(length(fits(c).vox),1)];
                
                
                [N,edges] = histcounts(y,boot.fitRange);
                [errorObj,lineObj,mn,se] = scatterline(x,y,boot.fitRange(N>boot.binThresh),NaN,1000,roiColors(ROI)*mult,2,1); hold on;
                extX = [0:.1:10];
                extLine = plot(extX, extendLine(lineObj,extX),'Color',roiColors(ROI)*mult,'LineWidth',.5);
                set(lineObj,'LineWidth',.5); set(extLine,'LineStyle',':');
                set(errorObj,'FaceAlpha',.2);
                
        end
        
        
        if c == 2
            xl = xlim; xlim([0 6]); y = ylim; if containsTxt(ROI,'faces') ylim([0 yl]); else ylim([0 6]);end
            set(gca,'TickDir','out');
            xlabel('Eccen (dva)','FontSize',fontSize); ylabel(['Size (SD/sqrt(N)) (dva)'],'FontSize',fontSize);
            axis square;
            title(ROI{1});
        end
        
        %superTitle(titleText,titleSize);
        

end


if saveFig == 1
    txt = ['crossVox_boot_' whichModel '_' whichStim '_r' num2str(minR2) '_' hemText(hems)];
    figure(f(1));
    niceSave([dirOf(pwd) 'figures/' expt '/sizeEccen/'],['scatter_' ROI{1} '_' txt '_' bootMethod '_dotAlpha' num2str(doAlpha)],[],[],{'png' 'svg'}); % just save pngs, since these can be generated pretty quickly

end
playSound;