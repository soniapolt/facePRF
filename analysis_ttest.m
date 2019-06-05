% loads & plots distributions of XY changes/anything else per each subject

clear all; close all;

subjs = {'TH' 'DF' 'EM' 'JG' 'MG' 'SP'};
expt = 'fixPRF';
tests = {'Y' 'X' 'gain' 'size'}; % can be parname (Y,X,sd,gain,exp,shift) or pRF.read value (r2,size,eccen,gain)
whichM = 'median'; % mean or median


minR2 = 20;          % cutoff for vox selection
ROIs= {'V1'};%{'hV4' 'IOG_faces' 'pFus_faces' 'mFus_faces'};%
fitSuffix = '';

whichStim = 'photo';%'eyes';%
whichModel = 'kayCSS';
hems = {'lh' 'rh'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems,fitSuffix));
ROInum = cellNum(ROIs,info.ROIs);
subjNum = cellNum(subjs,info.subjs);

comps = nchoosek(1:length(roi(1).fits),2);
fprintf('\n%s\n\n**************\n',whichM);
for t = 1:length(tests)
    test = tests{t};
    fprintf('-----\n%s:\n-----\n',test);
    for cc = 1:size(comps,1)
        baseCond = comps(cc,1);
        c = comps(cc,2);
        comp(c).descr = [roi(1).fits(baseCond).cond ' vs ' roi(1).fits(c).cond];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  ttest the parameters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        parNames = roi(1).fits(1).parNames;
        testNum = cellNum(test,parNames);
        
        
        
        for r = 1:length(ROIs)
            sData = [];
            for s = 1:length(subjs)
                % grab basecond
                vox = subj(subjNum(s)).roi(ROInum(r)).fits(baseCond).vox;
                % grab comp cond
                vox2 = subj(subjNum(s)).roi(ROInum(r)).fits(c).vox;
                
                if ~isempty(testNum)
                    pars = vertcat(vox.params);
                    sData(s,1) = nanmean(pars(:,testNum));
                    pars2 = vertcat(vox2.params);
                    sData(s,2) = nanmean(pars2(:,testNum));
                else
                    eval(['sData(s,1) = nan' whichM '([vox.' test ']);']);
                    eval(['sData(s,2) = nan' whichM '([vox2.' test ']);']);
                end
            end
            
            comp(c).groupData{r} = sData;
            [H,P,CI,STATS] = ttest(sData(:,1),sData(:,2));
            comp(c).diff = mean(sData(:,2))-mean(sData(:,1));
            if strcmp(test,'X') || strcmp(test,'Y')
                comp(c).diff = comp(c).diff /roi(1).fits(1).ppd; end
            comp(c).p(r) = P;
            comp(c).stats{r} = STATS;
            if H sig = '***'; else sig = ''; end
            fprintf('%s [%s %s:] %s param %s, %s: t(%d)=%.2f, p=%.3f; diff = %.3f\n',sig,hemText(hems),ROIs{r},test,comp(c).descr,whichM,STATS.df,STATS.tstat,P,comp(c).diff);
        end
    end
end
