% loads & plots distributions of XY changes/anything else per each subject

clear all; close all;

subjs = {'TH' 'DF' 'EM' 'JG' 'MG' 'SP'};
expt = 'fixPRF';
tests = {'Y' 'X' 'eccen' 'size' 'gain' 'r2'}; % can be parname (Y,X,sd,gain,exp,shift) or pRF.read value (r2,size,eccen,gain)
whichM = 'median'; % mean or median

minR2 = 50;          % cutoff for vox selection
ROIs= standardROIs('face');%{'V1' 'hV4' 'IOG_faces' 'pFus_faces' 'mFus_faces'};
fitSuffix = '';

whichStim = 'photo';
whichModel = 'kayCSS';
hems = {'lh' 'rh'};

factNames = {'hem' 'ROI' 'condition'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hem = struct;
for h = 1:length(hems)
    load(pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,{hems{h}},fitSuffix));
    hem(h).subj = subj;
end

ROInum = cellNum(ROIs,info.ROIs);
subjNum = cellNum(subjs,info.subjs);
fprintf('\n**************\n%s, %s\n**************\n',whichM, strTogether(ROIs));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ttest the parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parNames = roi(1).fits(1).parNames;

% data formatting for: RMAOV33.m
% each row = [data fact1# fact2# fact3# subj#];


for t = 1:length(tests)
    test = tests{t};
    anovaData = [];
    testNum = cellNum(test,parNames);
    
    for h = 1:length(hems)
        for r = 1:length(ROIs)
            sData = [];
            for s = 1:length(subjs)
                for c = 1:length(hem(1).subj(1).roi(1).fits)
                    % grab voxels
                    vox = hem(h).subj(subjNum(s)).roi(ROInum(r)).fits(c).vox;
                    
                    if ~isempty(testNum)
                        pars = vertcat(vox.params);
                        eval(['sData = nan' whichM '(pars(:,testNum));']);
                    else
                        eval(['sData = nan' whichM '([vox.' test ']);']);
                    end
                    % an = struct('factor',{'hem' 'ROI' 'condition'},'levels',{hems ROIs {'inverted' 'upright'}});
                    anovaData = [anovaData; sData h r c s];
                end
            end
        end
    end
    
    % RMAOV33(anovaData,.05,factNames);
    result = rmAnova3(anovaData,factNames,0);
    
    % quick ANOVA summary
    fprintf('\n-----\n%s:\n-----\n',test);
    main.s = []; main.ns = []; int.s = []; int.ns=[];
    for n = 1:length(result)
        
            text = sprintf('%s, F(%d,%d)=%.2f, p=%.3f.\n',result(n).name,result(n).df(1),result(n).df(2),result(n).F,result(n).p);
        if strcmp(result(n).type,'main')
            if result(n).h
        main.s = [main.s text];
        else
        main.ns = [main.ns text];end
        else
        if result(n).h
        int.s = [int.s text];else
       	int.ns = [int.ns text];end    
        end  
    end
    fprintf('Significant Main Effects: \n%s\n',main.s);
    fprintf('Significant Interactions: \n%s\n',int.s);
    fprintf('N.S: \n%s%s\n',main.ns,int.ns);
end
%             if H sig = '***'; else sig = ''; end
%             fprintf('%s [%s %s:] %s param %s, %s: t(%d)=%.2f, p=%.3f\n',sig,hemText(hems),ROIs{r},test,comp(c).descr,whichM,STATS.df,STATS.tstat,P);
%