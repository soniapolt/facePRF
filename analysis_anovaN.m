% loads & plots distributions of XY changes/anything else per each subject

clear all; close all;

subjs =prfSubjs;
expt = 'fixPRF';
tests = {'Y' 'X' 'eccen' 'size' 'gain' 'r2'}; % can be parname (Y,X,sd,gain,exp,shift) or pRF.read value (r2,size,eccen,gain)
whichM = 'median'; % mean or median

minR2 = 20;          % cutoff for vox selection
ROIs= standardROIs('face-');%{'V1' 'hV4' 'IOG_faces' 'pFus_faces' 'mFus_faces'};
fitSuffix = '';

txtName = 'anovan';

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
checkDir('results/');
    fid = fopen(['results/ANOVA_' txtName '_' whichModel '_' whichStim],'w+');

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
                    if ~isempty(vox)
                    if ~isempty(testNum)
                        pars = vertcat(vox.params);
                        eval(['sData = nan' whichM '(pars(:,testNum));']);
                    else
                        eval(['sData = nan' whichM '([vox.' test ']);']);
                    end
                    
                    % an = struct('factor',{'hem' 'ROI' 'condition'},'levels',{hems ROIs {'inverted' 'upright'}});
                    anovaData = [anovaData; sData h r c s];
                    else error('Missing data! Can''t run this ANOVA function!'); end
                end
            end
        end
    end
    
    [p,t,stats,terms] = anovan(anovaData(:,1),anovaData(:,2:end-1),'varnames',factNames,'display','off');
    
    % quick ANOVA summary
    
    fprintf('\n-----\n%s:\n-----\n',test);
    fprintf(fid,'\n-----\n%s:\n-----\n',test);
    
    main.s = []; main.ns = []; int.s = []; int.ns=[];
    for n = 1:length(t)
        
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
    fprintf(fid,'Significant Main Effects: \n%s\n',main.s);
    fprintf(fid,'Significant Interactions: \n%s\n',int.s);
    fprintf(fid,'N.S: \n%s%s\n',main.ns,int.ns);
    
    fprintf('Significant Main Effects: \n%s\n',main.s);
    fprintf('Significant Interactions: \n%s\n',int.s);
    fprintf('N.S: \n%s%s\n',main.ns,int.ns);
end
fclose(fid);