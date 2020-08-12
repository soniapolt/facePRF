% loads & plots distributions of XY changes/anything else per each subject
% not used in manuscript figures

clear all; close all;

subjs = prfSubjs;%{'MG' 'JG' 'TH' 'EM' 'DF' 'SP'};%
expt = 'fixPRF';

saveFig =1;
convertDVA = 1;
flipX = 1;

minR2 = 'r2-20';          % cutoff for vox selection
ROIs= ['V1' standardROIs('face')];%{'mFus_faces'};%

whichStim = 'outline';%'photo';%'internal';%'eyes';%
whichModel = 'kayCSS';%'inflipCSSn';%'flipCSSn';%'cssExpN';%'cssShift';%
whichM = 3; % 1 = mean, 2 = mode/peak, 3 = median
mS = {'mean' 'mode' 'median'};
plotPars = {'Y'};%{'size' 'X' 'Y' 'gain'};%{'gain' 'r2' 'Y' 'X' };
parTitles = {'Y Estim'};%{'Size [Sigma/sqrt(N)] (dva)' 'X Estim.' 'Y Estim.' 'Gain Estim'};%{'Gain Estim' 'Estimated R^{2}' 'Y Estim.' 'X Estim.' 'Size [Sigma/sqrt(N)] (dva)'};
plotType = {'colorbar'};%{'box' 'box' 'box' 'box'};%{'box' 'distr' 'distr' 'distr' 'box'}; % bars, distr, or delta

hems = {'lh' 'rh'};
fitSuffix = '';%'_orig';%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 11; titleSize = 14;
nBins = 20; % histogram bins

% now we load in the data from both hemispheres, and threshold across
load(pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems,fitSuffix));
ROInum = cellNum(ROIs,info.ROIs);
subjNum = cellNum(subjs,info.subjs);

if length(subjNum) == 1 roi = subj(subjNum).roi; else subj = subj(subjNum); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1: distribution of parameters for this ROI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p = 1:length(plotPars)
    
    titleText = [whichModel ' ' parTitles{p} ', Subj: '];
    titleText = [titleText strTogether(subjs) ' (voxels R^2 > ' num2str(minR2) '), ' whichStim ' stim, ' whichModel ' model'];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if saveFig && onLaptop figSize = [.1 .1 .8 .8]; else figSize = [.2 .1 8 .8]; end
    niceFig(figSize,fontSize);
    numPlots = [2 ceil(length(ROIs)/2)];pl = 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for r = 1:length(ROIs)
        
        

        sPars = []; mPars = [];
        for c = 1:length(roi(1).fits)
            for s = 1:length(subjNum)
                fits = subj(subjNum(s)).roi(ROInum(r)).fits(c);
                try
                cPars{c} = getPar(plotPars{p},fits,1,flipX);
                catch cPars{c} = NaN; end % missing ROIs
                
                % aggregate across subjects
                % correction for size def
                if containsTxt(plotPars{p},'size')
                cPars{c}= cPars{c}./2;    
                end
                sPars{s,c} = cPars{c}; % full distribution for this subject
                
                eval(['mPars(s,c) = nan' mS{whichM} '(cPars{c});']); % mean value for this subject
            end
        end
        hues = [1 .5];
        
        switch plotType{p}
            case 'bars'
                niceBars2(mPars,mS{whichM},1,{fits.cond},[roiColors(ROIs(r)).*hues(1); roiColors(ROIs(r)).*hues(2)]);
            case 'delta'
                niceBars2(mPars(:,2)-mPars(:,1),mS{whichM},1,[{fits(2).cond ' - ' fits(1).cond}],roiColors(ROIs(r)));
            case 'distr'
                if containsTxt(plotPars{p},'size')
                    xl = [0 10] ;
                else xl = [-5 5]; end
                %plotMeanDistr(cData,nBins,color,whichM)
                
                for c = 1:length(roi(1).fits)
                    [h,meds, normcounts] = plotMeanDistr(sPars(:,c),nBins,roiColors(ROIs(r)).*hues(c),1);xlim(xl);
                end
                %ylim([0 nanmax(normcounts(:))]);
                %                 sample1 = sPars{:,1};
                %                 [pv, od, effectsize] = permutationTest(sPars{:,1}, sPars{:,2}, 1000), ...
                %                      'plotresult', 1, 'showprogress', 0);
                %                  fprintf('%s Permutation Test, %s: p = %.6f, observed difference = %.2f\n ',...
                %                      whichModel, plotPars{p},pv,od);
            case 'deltadist'
                if containsTxt(plotPars{p},'size')
                    xl = [0 10] ;
                else xl = [-5 5]; end
                for n = 1:length(sPars)
                    dPars{n} = sPars{n,2}-sPars{n,1};
                end
                [h,peak, normcounts] = plotMeanDistr(dPars,nBins,roiColors(ROIs(r))); xlim(xl);
            case 'box'
                if containsTxt(plotPars{p},'gain')
                    yl = [0 5];
                elseif containsTxt(plotPars{p},'r2')
                    yl = [0 100];
                elseif containsTxt(plotPars{p},'size')
                    yl= ([0 5]);
                else yl = [];
                end
                niceBoxplot(fliplr(mPars),fliplr({roi(1).fits(1).cond roi(1).fits(2).cond}),1,flipud([roiColors(ROIs(r)).*hues(1); roiColors(ROIs(r)).*hues(2)]),yl);
            case 'colorbar'
                
                
            
        end
        set(gca,'TickDir','out');
        pl = pl+1;
        axis square;
        
        [h,pv,ci] = ttest(mPars(:,2)-mPars(:,1));
        title({ROIs{r};parTitles{p};['ttest on medians: p=' num2str(pv)]},'fontSize',titleSize,'interpreter','none','FontWeight','bold');
        %xlabel(roi(1).fits(1).parNames{p},'fontSize',titleSize);
    end
    superTitle(titleText,titleSize,.025);
    
    if saveFig == 1
        if length(subjs) == 1 txt = ['subj' subjs{1}]; else txt = ['groupN' num2str(length(subjs))]; end
        txt = [plotPars{p} '_' hemText(hems) '_' txt '_flipX' num2str(flipX)];
        niceSave([dirOf(pwd) 'figures/' expt '/crossSubj/'],txt,[],[],{'png' 'svg'}); % just save pngs, since these can be generated pretty quickly
    end
end
if onLaptop playSound; end