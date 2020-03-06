% loads & plots some basics re: invPRF fits of the cssModel (positive pRF
% constraint) & the cssShift model (added a baseline term)
%
% SP 4/25/18
% this version plots ROIs as rows, params as columns, subjects as separate graphs
%
% version 2 implements standard pRF trimming, which we can hopefully use
% across different visualizations of the fits

clear all; close all;

subjs = prfSubjs;%{'DF' 'EM' 'TH' 'JG' 'MG' 'SP'};%;%;
task = '';
expt = 'fixPRF';
saveFig = 1;
convertDVA = 1;

whichModels = {'kayCSS'};% 'kayCSS','cssExpN' 'cssExpN' 'intempCSSn' 'inflipCSSn'};
modelStims = {'outline'};%,'internal', 'outline','internal','internal','internal'};

fitSuffix = '';%'_orig';%

fitsSuffix = ''; %'_orig';
compConds = [2 1];

minR2 = ['r2-20'];
ROI=standardROIs(7);
whichM = 'median';


plotPars = {'Y'};%{'r2' 'Y' 'size' 'gain'};
parTitles = {'Estimated R^{2}' 'Y Estim.' 'Size [2xSD/sqrt(N)] (dva)' 'Gain Estim'};
%plotType = {'fit' 'fit' 'fit' 'fit'}; % 1 = boxplot, 2 = distr, 3 = scatter, 4 = delta(distribution), 5 = histfit (delta)
%plotType = {'scatter' 'distr' 'distr' 'box'}; % 1 = boxplot, 2 = distr, 3 = scatter, 4 = delta(distribution)
plotType = {'delta' 'delta' 'delta' 'delta'}; % 1 = boxplot, 2 = distr, 3 = scatter, 4 = delta(distribution)

hems = {'rh' 'lh'};
fitSuffix = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fontSize = 14; titleSize = 18;




for p = 1:length(plotPars)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if saveFig && onLaptop figSize = [0 0 1 1]; else figSize = [.1 .1 .9 .9]; end
    f(p) = niceFig(figSize,fontSize);
    [numPlots(1) numPlots(2)] = subplotDims(length(modelStims),[],2); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

for t = 1:length(modelStims)
    
    % load this model
    whichModel = whichModels{t};
    load(pRFfile(dirOf(pwd),expt,minR2,modelStims{t},whichModel,hems,fitSuffix,task));
    ROInum = cellNum(ROI,info.ROIs);
    subjNum = cellNum(subjs,info.subjs);
    if length(subjNum) == 1 roi = subj(subjNum).roi; end
    % assuming either one subject, or all subjects
    fits = roi(ROInum).fits;
    
    for p = 1:length(plotPars)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(f(p)); subplot(numPlots(1),numPlots(2),t)
        % get param values for this condition
        
        for cc = 1:length(compConds)
            c = compConds(cc);
            parNum = cellNum(plotPars{p},fits(1).parNames);
            plPars{cc} = getPar(plotPars{p},fits(c))
%             if ~isempty(parNum)
%                 pars = vertcat(fits(c).vox.params);
%                 plPars{cc} = pars(:,parNum)';
%             else
%                 eval(['plPars{cc} = [fits(c).vox.' plotPars{p} '];']);  end
%             
%             if convertDVA && containsTxt(plotPars{p},'Y') || containsTxt(plotPars{p},'X') || containsTxt(plotPars{p},'sd')
%                 % rescale some parameters so that they are in DVA units and
%                 % centered around zero (center of screen)
%                 if ~containsTxt(plotPars{p},'sd') % don't re-center the SD
%                     plPars{cc} = fits(1).res-plPars{cc}-roi(1).fits(1).res/2;
%                 end
%                 plPars{cc} = plPars{cc}./roi(1).fits(1).ppd;
%             end
        end
        
        switch plotType{p}
            case 'box'
                if containsTxt(plotPars{p},'gain')
                    niceBoxplot([plPars{1};plPars{2}]',{fits(compConds(1)).cond fits(compConds(2)).cond},1,[condColors(4);condColors(2)],[0 10]);
                    ylim([0 5]);
                elseif containsTxt(plotPars{p},'r2')
                    niceBoxplot([plPars{1};plPars{2}]',{fits(compConds(1)).cond fits(compConds(2)).cond},1,[condColors(4);condColors(2)],[minR2 100]);
                else
                    niceBoxplot([plPars{1};plPars{2}]',{fits(compConds(1)).cond fits(compConds(2)).cond},1,[condColors(4);condColors(2)]);
                end
                if containsTxt(plotPars{p},'size')
                    ylim([0 5]);
                end
            case 'distr'
                nBins = 20;
                plotDistr(plPars,1,{fits(compConds(1)).cond fits(compConds(2)).cond},nBins,3,1);
                if containsTxt(plotPars{p},'size')
                    xlim([0 10]);
                else xlim([-5 5]); end
                
                [pv, od, effectsize] = permutationTest(plPars{1}, plPars{2}, 1000, ...
                    'plotresult', 0, 'showprogress', 0);
                fprintf('%s Permutation Test, %s: p = %.6f, observed difference = %.2f\n ',...
                    whichModel, plotPars{p},pv,od);
            case 'delta'
                 %v = vline(0,'k'); set(v,'LineWidth',2); hold on;
                 niceHist(plPars{2}-plPars{1},condColors(c),1);
                 [~,pval,~,stats] = ttest(plPars{2}-plPars{1});
                 xlabel({[plotPars{p} ' delta' ];sprintf('t(%d) = %0.2f, p = %0.3f',stats.df,stats.tstat,pval)},'FontSize',10);
                 ylabel('Count');
            case 'fit'
                numBins = 100; plotMu = 1;
                [mu, ci] =  niceHistFit(plPars{2}-plPars{1},condColors(c),numBins,plotMu);
                xlabel('Delta'); ylabel('Count');
                xlabel({[plotPars{p} ' delta' ];sprintf('mu = %0.2f (CI = [%0.2f%0.2f])',mu,ci(1),ci(2))},'FontSize',10);
                
                
            case 'scatter'
                hold on;
                scatterCent(plPars{1},plPars{2},condColors(c),...
                    fits(compConds(1)).cond,fits(compConds(2)).cond,[],fontSize,0,1);
                if strcmp(plotPars{p},'r2') xlim([0 100]); ylim([0 100]); end
                hold on; xl = get(gca,'xlim');
                hold on; plot(xl,xl,'k:');
        end
        title([whichModel ' ' modelStims{t}],'FontSize',titleSize);
        axis square;
    end
    
end

for p = 1:length(plotPars)
figure(f(p));
superTitle([ROI{1} ' ' upper(plotPars{p})],titleSize,.97);
superTitle(sprintf('%s Model Comparisons, %s %s, R2 cutoff = %s',ROI{1},whichM,plotPars{p},minR2),titleSize-4,.05);

%subplotresize;
if saveFig == 1
    if length(subjs) == 1
        txt = ['subj' subjs{1}];
    else txt = ['groupN' num2str(length(subjs))]; end
    txt = [hemText(hems) '_' ROI{1} '_modelComp_' txt '_' plotPars{p}];
    if isequal(plotType{1:end}) txt = [txt '_' plotType{1}]; end
    niceSave([dirOf(pwd) 'figures/' expt '/modelComp/'],txt); % just save pngs, since these can be generated pretty quickly
end
end
if onLaptop playSound; end