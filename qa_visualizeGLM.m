clear all; %close all;

subj = 'JP';
sessNum = 1;
expt = 'fixPRF';

switch expt
            case 'invPRF3'
                [session, numRuns] = vpnlSessions(expt,subj,sessNum,'fix');
                numTRs = 186;
            case 'fixPRF'
                [session, numRuns] = vpnlSessions(expt,subj);
                numTRs = 136;
            case 'compPRF'
                [session, numRuns] = vpnlSessions(expt,subj,sessNum);
                numTRs = 186;
end
        %numRuns =1;
%if containsText(subj,'SP') session = [session '_0']; end

whichStim = 'photo';
whichModel = 'kayCSS';

ROI ='rh_mSTS_faces';%'rh_V1';%'%%'lh_V2';%

voxNum =1; % if 0, random
sp1=4; sp2=ceil(numRuns/sp1);


[dataName, fitsName] = fitsDirs(dirOf(pwd),expt,session,whichStim,vpnlROI(ROI,subj),whichModel);
load(dataName);
fprintf('Loading data: %s...\n',dataName);

[val,sortVox] = sort(mv.glm.varianceExplained,2,'descend');


if voxNum==0
    voxNum = randi(size(mv.tSeries,2)); 
else voxNum = sortVox(voxNum); end

% sortVox(isnan(val))
% voxNum = sortVox(3);
TCs = mv.tSeries(:,voxNum);
TCs = reshape(TCs,numRuns,numTRs);

resids = mv.glm.residual(:,voxNum);
resids = reshape(resids,numRuns,numTRs);

niceFig([.1 .1 .8 .8],14);
for n = 1:numRuns
subplot(sp1,sp2,n)
plot(TCs(n,:)); xlim([1 numTRs]);
hold on; plot(resids(n,:),'r');
legend('TC','resids');
end
superTitle(['Vox: ' num2str(voxNum) ', Var Explained: ' num2str(mv.glm.varianceExplained(voxNum))]);

DMs = reshape(mv.glm.designMatrix,[numTRs,numRuns,size(mv.glm.designMatrix,2)]);

niceFig([.1 .1 .6 .5],14);
for n = 1:numRuns
subplot(sp1,sp2,n)
imshow(squeeze(DMs(:,n,:))')
end
superTitle('Design Matrices');

%imshow(mv.glm.designMatrix')