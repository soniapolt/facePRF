

x = [-50:.1:50];

%%% CSS model


% the parameters of the CSS model are [R C S G N] where
%   R is the row index of the center of the 2D Gaussian
%   C is the column index of the center of the 2D Gaussian
%   S is the standard deviation of the 2D Gaussian
%   G is a gain parameter
%   N is the exponent of the power-law nonlinearity

%modelfun = @(pp,dd) pp(4)*((dd*vflatten(makegaussian2d(stimSize,pp(1),pp(2),pp(3),pp(3),xx,yy,0,0)/(2*pi*pp(3)^2))).^pp(5));


%sigmas = [1]; exps = [.2 1]; gain = 1; 

% l1norm is used to make the sum of the gaussian == 1

%for s = 1:length(sigmas), sigma = sigmas(s)
   % for e = 1:length(exps), exp = exps(e)
   
   
sigma = 1; gain = 1;

% line 1 - nonCSS model
exp = 1;
l1norm = (2*pi*sigma^2);
prf = gain*((normpdf(x,0,sigma)/l1norm).^exp);
sz = sigma*sqrt(exp);
plot(x,prf,'r'); hold on; 
vline(-sz,'r:'); vline(sz,'r:',num2str(sz));

% line 2 - css model
exp = .2;
l1norm = (2*pi*sigma^2);
prf = gain*((normpdf(x,0,sigma)/l1norm).^exp);
sz = sigma*sqrt(exp);
plot(x,prf,'g'); hold on; 
vline(-sz,'g:'); vline(sz,'g:',num2str(sz));

% line 3 - inserting css size as sigma (like i maybe did incorrectly in the
% optimRec simulation)
exp = .2; sigma = sigma/sqrt(exp);
l1norm = (2*pi*sigma^2);
prf = gain*((normpdf(x,0,sigma)/l1norm).^exp);
sz = sigma*sqrt(exp);
plot(x,prf,'b'); hold on; 
vline(-sz,'b:'); vline(sz,'b:',num2str(sz));
    %end
%end

%%% this is what i used for the optimRec simulation...
[XX,YY]=meshgrid(-50:50);
cov= PRF(XX,YY,0,0,fits(c).vox(v).size*sim.ppd);
