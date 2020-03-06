% loads & plots distributions of XY changes/anything else per each subject
% subscript to look at effects of condition and hemisphere within ROIs

clear all; close all;

subjs =prfSubjs;
expt = 'fixPRF';
tests = {'Ydeg' 'Xdeg' 'eccen' 'size' 'gain' 'r2'}; % can be parname (Y,X,sd,gain,exp,shift) or pRF.read value (r2,size,eccen,gain)
whichM = 'median'; % mean or median

r2cutoff = 'r2-20';%'perc-50';%;%  %    %or 'r2-20'    % cutoff for vox selection

ROIs= standardROIs;%{'V1' 'hV4' 'IOG_faces' 'pFus_faces' 'mFus_faces'};
fitSuffix = '';


whichStim = 'outline';%'edge';%'photo';%'internal';%
whichModel = 'kayCSS';%'intempCSSn';%
hems = {'lh' 'rh'};

factNames = {'condition' 'hem'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hem = struct;

for h = 1:length(hems)
    pF = pRFfile(dirOf(pwd),expt,r2cutoff,whichStim,whichModel,{hems{h}},fitSuffix);
    load(pF); fprintf('Loading %s\n...',pF);
    hem(h).subj = subj;
end

ROInum = cellNum(ROIs,info.ROIs);
subjNum = cellNum(subjs,info.subjs);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ttest the parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parNames = roi(1).fits(1).parNames; checkDir([dirOf(pwd) '/stats/' expt '/intraROI']);

for r  = 1:length(ROIs)
% data formatting for: RMAOV33.m
% each row = [data fact1# fact2# fact3# subj#];
    outputFile = [dirOf(pwd) 'stats/' expt '/intraROI/ANOVA2_' ROIs{r} '_' whichModel '_' whichStim '_' r2cutoff '.txt'];
    fid = fopen(outputFile,'w+');
fprintf(fid,'\n**************\n%s, %s\n**************\n',whichM, ROIs{r});
fprintf(fid,['pRF file: %s\n'],pF,ROIs{r});
fprintf(fid,'\n**************\n%s, %s\n**************\n',whichM, ROIs{r});


for t = 1:length(tests)
    test = tests{t};
    anovaData = [];
    testNum = cellNum(test,parNames);
    
    rmSubjs = [];
            sData = [];
            for s = 1:length(subjs)
                for c = 1:length(subj(1).roi(1).fits)
                    for h = 1:length(hems)
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
                    anovaData = [anovaData; sData c h subjNum(s)];
                    %else error('Missing data! Can''t run this ANOVA function!'); end
                    else rmSubjs(end+1) =  subjNum(s); if c == 1 
                            fprintf('Missing data from %s in %s %s! \n', subjs{s}, hems{h}, ROIs{r});
                            fprintf(fid,'Missing data from %s in %s %s!\n', subjs{s}, hems{h}, ROIs{r});
                        end
                    end
                end
                end
            end
    
    % check for missing data and remove those subjects from the comparison
    
    if ~isempty(rmSubjs) for s = unique(rmSubjs);
        anovaData(find(anovaData(:,end)==s),:) = []; 
    fprintf('Removed data from subj %s...\n',subjs{subjNum(s)}); end
    end
    
    % % function stats = rm_anova2(Y,S,F1,F2,FACTNAMES)
    result = rm_anova2(anovaData(:,1),anovaData(:,end),anovaData(:,2),anovaData(:,3),factNames);
    
    % quick ANOVA summary
    
    fprintf('\n-----\n%s %s:\n-----\n',ROIs{r},test);
    fprintf(fid,'\n-----\n%s %s:\n-----\n',ROIs{r},test);
    
    main.s = []; main.ns = []; int.s = []; int.ns=[];
    
    for n = 2:4
        
            text = sprintf('%s, F(%d)=%.2f, p=%.3f.\n',result{n,1},result{n,3},result{n,5},result{n,6});
        if ~containsTxt(result{n,1},' x ')
            if result{n,6} < .05
        main.s = [main.s text];
        else
        main.ns = [main.ns text];end
        else
        if result{n,6} < .05
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
fprintf('Saving %s.\n',outputFile);
end
if onLaptop playSound; end
