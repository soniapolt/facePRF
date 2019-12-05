% loads & plots distributions of XY changes/anything else per each subject

clear all; close all;

subjs = prfSubjs;%{'MG' 'JG' 'TH' 'EM' 'DF' 'SP'};%
expt = 'fixPRF';

saveFig =1;
convertDVA = 1;

minR2 = 'r2-20';          % cutoff for vox selection
ROIs= standardROIs('face+');%;%['hV4' standardROIs('face+')]

whichStim = 'outline';%'photo';%'internal';%'eyes';%
whichModel = 'kayCSS';%'inflipCSSn';%'flipCSSn';%'cssExpN';%'cssShift';%
whichM = 3; % 1 = mean, 2 = mode/peak, 3 = median
mS = {'mean' 'mode' 'median'};
plotPars = {'gain' 'r2' 'Y' 'X' 'size'};{'size'};%
parTitles = {'Gain Estim' 'Estimated R^{2}' 'Y Estim.' 'X Estim.' 'Size [2xSD/sqrt(N)] (dva)'};{'Size Estim.'};%
plotType = {'box' 'distr' 'distr' 'distr' 'box'}; % bars, distr, or delta

hems = {'rh' 'lh'};
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
    if saveFig && onLaptop figSize = [0 0 1 1]; else figSize = [.2 .1 8 .8]; end
    niceFig(figSize,fontSize);
    numPlots = [2 ceil(length(ROIs)/2)];pl = 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for r = 1:length(ROIs)
        
        subplot(numPlots(1),numPlots(2),pl)
        
        sPars = []; mPars = [];
        for c = 1:length(roi(1).fits)
            for s = 1:length(subjNum)
                fits = subj(subjNum(s)).roi(ROInum(r)).fits;
                
                parNum = cellNum(plotPars{p},fits(1).parNames);
                if ~isempty(parNum)
                    pars = vertcat(fits(c).vox.params);
                    if ~isempty(pars)
                        cPars{c} = pars(:,parNum)';
                    else cPars{c} = NaN;
                    end
                else
                    eval(['cPars{c} = [fits(c).vox.' plotPars{p} '];']);
                    if isempty(cPars{c}) cPars{c} = NaN; end
                end
                
                
                % rescale some parameters so that they are in DVA units and
                % centered around zero (center of screen)
                if convertDVA && containsTxt(plotPars{p},'Y') || containsTxt(plotPars{p},'X') || containsTxt(plotPars{p},'sd')
                    if ~containsTxt(plotPars{p},'sd') % don't re-center the SD
                        cPars{c} = fits(1).res-cPars{c}-roi(1).fits(1).res/2;
                    end
                    cPars{c} = cPars{c}./roi(1).fits(1).ppd;
                end
                
                % aggregate across subjects
                sPars{s,c} = cPars{c}; % full distribution for this subject
                eval(['mPars(s,c) = nan' mS{whichM} '(cPars{c});']); % mean value for this subject
            end
        end
        hues = [1 .5];
        switch plotType{p}
            case 'bars'
                niceBars2(mPars,mS{whichM},1,{fits.cond},[condColors(r,1).*hues(1); condColors(r,1).*hues(2)]);
            case 'delta'
                niceBars2(mPars(:,2)-mPars(:,1),mS{whichM},1,[{fits(2).cond ' - ' fits(1).cond}],condColors(r,1));
            case 'distr'
                if containsTxt(plotPars{p},'size')
                    xl = [0 10] ;
                else xl = [-5 5]; end
                %plotMeanDistr(cData,nBins,color,whichM)
                
                for c = 1:length(roi(1).fits)
                    [h,peak, normcounts] = plotMeanDistr(sPars(:,c),nBins,condColors(r,1).*hues(c));
                end
                %ylim([0 nanmax(normcounts(:))]);
                %                 sample1 = sPars{:,1};
                %                 [pv, od, effectsize] = permutationTest(sPars{:,1}, sPars{:,2}, 1000), ...
                %                      'plotresult', 1, 'showprogress', 0);
                %                  fprintf('%s Permutation Test, %s: p = %.6f, observed difference = %.2f\n ',...
                %                      whichModel, plotPars{p},pv,od);
            case 'box'
                if containsTxt(plotPars{p},'gain')
                    yl = [0 5];
                elseif containsTxt(plotPars{p},'r2')
                    yl = [0 100];
                elseif containsTxt(plotPars{p},'size')
                    yl= ([0 5]);
                else yl = [];
                end
                niceBoxplot(mPars,{fits(1).cond fits(2).cond},1,[condColors(r,1).*hues(1); condColors(r,1).*hues(2)],yl);
        end
        pl = pl+1;
        axis square;
        
        [h,pv,ci] = ttest(mPars(:,2)-mPars(:,1));
        title({ROIs{r};parTitles{p};['ttest on medians: p=' num2str(pv)]},'fontSize',titleSize,'interpreter','none','FontWeight','bold');
        %xlabel(roi(1).fits(1).parNames{p},'fontSize',titleSize);
    end
    superTitle(titleText,titleSize,.025);
    
    if saveFig == 1
        if length(subjs) == 1 txt = ['subj' subjs{1}]; else txt = ['groupN' num2str(length(subjs))]; end
        txt = [plotPars{p} '_' hemText(hems) '_' txt];
        niceSave([dirOf(pwd) 'figures/' expt '/randomEffects/'],txt); % just save pngs, since these can be generated pretty quickly
    end
end
if onLaptop playSound; end