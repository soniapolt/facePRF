% loads & plots distributions of XY changes/anything else per each subject

clear all; close all;

subjs =prfSubjs;
expt = 'fixPRF';
tests = {'Y' 'X' 'eccen' 'size' 'gain' 'r2'}; % can be parname (Y,X,sd,gain,exp,shift) or pRF.read value (r2,size,eccen,gain)
whichM = 'median'; % mean or median

r2cutoff = 'r2-20'%;'perc-33';%%  %    %or 'r2-20'    % cutoff for vox selection
whichTest = 'face';

ROIs= standardROIs(whichTest);%{'V1' 'hV4' 'IOG_faces' 'pFus_faces' 'mFus_faces'};
fitSuffix = '';

txtName = [whichTest '-' r2cutoff];

whichStim = 'outline';%'edge';%'photo';%'internal';%
whichModel = 'kayCSS';%'intempCSSn';%
hems = {'lh' 'rh'};

factNames = {'ROI' 'condition'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pF = pRFfile(dirOf(pwd),expt,r2cutoff,whichStim,whichModel,hems,fitSuffix);
    load(pF); fprintf('%s\n...',pF);

ROInum = cellNum(ROIs,info.ROIs);
subjNum = cellNum(subjs,info.subjs);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ttest the parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parNames = roi(1).fits(1).parNames;

% data formatting for: RMAOV33.m
% each row = [data fact1# fact2# fact3# subj#];
checkDir([pwd '/results']);
    fid = fopen([pwd '/results/ANOVA2_' txtName '_' whichModel '_' whichStim],'w+');
fprintf(fid,'\n**************\n%s, %s\n**************\n',whichM, strTogether(ROIs));
fprintf(fid,['pRF file: %s\nROIs: ' repmat('%s ',1,length(ROIs)) '\n'],pF,ROIs{:});

fprintf('\n**************\n%s, %s\n**************\n',whichM, strTogether(ROIs));


for t = 1:length(tests)
    test = tests{t};
    anovaData = [];
    testNum = cellNum(test,parNames);
    
    rmSubjs = [];
        for r = 1:length(ROIs)
            sData = [];
            for s = 1:length(subjs)
                for c = 1:length(subj(1).roi(1).fits)
                    % grab voxels
                    vox = subj(subjNum(s)).roi(ROInum(r)).fits(c).vox;
                    if ~isempty(vox)
                   
                     if ~isempty(testNum)
                        pars = vertcat(vox.params);
                        eval(['sData = nan' whichM '(pars(:,testNum));']);
                    else
                        eval(['sData = nan' whichM '([vox.' test ']);']);
                    end
                    
                    % an = struct('factor',{'hem' 'ROI' 'condition'},'levels',{hems ROIs {'inverted' 'upright'}});
                    anovaData = [anovaData; sData r c subjNum(s)];
                    %else error('Missing data! Can''t run this ANOVA function!'); end
                    else rmSubjs(end+1) =  subjNum(s); if c == 1 fprintf('Missing data in %s %s-%s!\n', subjs{s}, ROIs{r});
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
    result = rm_anova2(anovaData(:,1),anovaData(:,end),anovaData(:,2),anovaData(:,3),factNames)
    
    % quick ANOVA summary
    
    fprintf('\n-----\n%s:\n-----\n',test);
    fprintf(fid,'\n-----\n%s:\n-----\n',test);
    
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