% loads & plots distributions of XY changes/anything else per each subject

clear all; close all;

subjs = {'TH' 'DF' 'EM' 'JG' 'MG' 'SP'};
expt = 'fixPRF';
test = 'Y'; % 'X' 'eccen' 'size' 'gain' 'r2' can be parname (Y,X,sd,gain,exp,shift) or pRF.read value (r2,size,eccen,gain)

minR2 = 20;          % cutoff for vox selection
ROIs= standardROIs;%{'V1' 'hV4' 'IOG_faces' 'pFus_faces' 'mFus_faces'};
fitSuffix = '';

whichStim = 'photo';
whichModel = 'kayCSS';
hems = {'lh' 'rh'};

%factNames = {'hem' 'ROI' 'condition'};
predNames = {'condition' 'r2' 'gain'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems,fitSuffix));

ROInum = cellNum(ROIs,info.ROIs);
subjNum = cellNum(subjs,info.subjs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ttest the parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parNames = roi(1).fits(1).parNames;

% data formatting for stepwiselm: stepwiselm(X,Y), where X is column
% vector Y as a response variable and the columns of the matrix X as
% predictor variables

    testNum = cellNum(test,parNames);
    for r = 1:length(ROIs)
        fprintf('\n**************\n%s stepwise regression, %s:\n**************\n', ROIs{r},test);

        stepData = [];
        voxNum = length(roi(ROInum(r)).fits(1).vox);
        stepPred = zeros(voxNum*length(roi(r).fits),length(predNames));
        
        for c = 1:length(roi(1).fits)
            % grab voxels
            vox = roi(ROInum(r)).fits(c).vox;
            pars = vertcat(vox.params);
            
            if ~isempty(testNum)
                sData = pars(:,testNum)';
            else
                eval(['sData = [vox.' test '];']);
            end
            stepData = [stepData sData];
            
            % parse predictors
            catVars = [];
            for n = 1:length(predNames)
                if strcmp(predNames{n},'condition')
                    pData = ones(1,voxNum)*c;
                    catVars = [catVars n];
                else
                    predNum = cellNum(predNames{n},parNames);
                    if ~isempty(predNum)
                        pData = pars(:,predNum);
                    else
                        eval(['pData = [vox.' predNames{n} '];']);
                    end
                    
                end
                ind = voxNum*(c-1)+1:voxNum*(c-1)+voxNum;
                stepPred(ind,n) = pData;
            end
        end
        
        varNames = predNames; varNames{end+1} = test;
        lm = stepwiselm(stepPred,stepData,'VarNames',varNames,'CategoricalVars',catVars,'Verbose',2)
 
    end


%     %