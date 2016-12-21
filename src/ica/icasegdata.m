% function icadata = icasegdata(datain,varargin)
%
% This function performs an ICA decomposition of segmented data. It is
% reccomended, but not required, to first screen for unique artifacts
% or amplifier artifacts using artscreen.m.  The more homogeneous the
% input data, the better the separation that can be achieved.
%
% The returned data structure has all of the original fields, plus now
% includes mixing and unmixing matrices, as well as amplitude spectra
% for each component and the percentage of variance that each component
% accounts for.  The power spectra are weighted by the trial variance
% so that very sparsely active artifactual components do not have
% 'normal-looking' spectra due to averaging.
%
% Required Input:
%    Data structure containing the fields .data and .sr 
%
% Optional Arguments: 
%
%    ncomps-  the number of components to solve for. Must be <= the number 
%             of channels and have > ncomps^3 good samples of
%             data. Default: min([60 nchans round(goodsamps^(1/3))])
%
%    nkeep-   the number of components to retain for review. The comps
%             are sorted by variance, thus only minor artifacts comprise
%             the later components. Defaut: min([40 ncomps])
%
%    fftfreq- the max frequency to calculate for the component amplitude
%             spectra.  Default: 50
%
%    extica-  flag to use the 'extended' infomax ICA, which can separate
%             sub-gaussian sources such as line noise. Note that the
%             decomposition will take longer. Default: 1
%
%    badchans - channels that are not zeroed out but should be considered
%               bad and excluded from the ICA.
%                    
%
% Part 2 of artscreenEEG's basic data cleaning functions:
%    artscreen.m => icasegdata.m => icareview.m => icatochan.m (optional)
%
% Copyright (C) 2013 Cort Horton, <chorton@uci.edu>
% Copyright (C) 2016 Michael D. Nunez, <mdnunez1@uci.edu>
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

function datain = icasegdata(datain,varargin)

% VERSION HISTORY:
%   1.0 - created by Cort Horton  6/26/12
%   1.2 - Some changes to default values 7/1/13
%   1.3 - Now calculates the percentage of the variance in the data 
%          that each component accounts for.  7/15/13    
%   1.4 - Changed the default number of components to solve for. 8/23/16 -Michael Nunez 
%   1.5 - Export citation to the terminal when Infomax ICA. 12/20/16 - Michael Nunez

%To do:
% 1) Find fastica alogirthm, add fastica as an option
% 2) Add SOBI as an option
% 3) Change default to KEEP random components?

if nargin < 1; help icasegdata; return; end;

% Parse inputs;
[~,ncomps,nkeep,fftfreq,extica,badchans]=...
    parsevar(varargin,'ncomps',[],'nkeep',[],...
    'fftfreq',50,'extica',1,'badchans',[]);

fprintf('Infomax ICA used! Please cite:\n');
fprintf('\n');
fprintf('Makeig, S., Bell, A.J., Jung, T-P and Sejnowski, T.J., \n');
fprintf('Independent component analysis of electroencephalographic data.\n');
fprintf('In: D. Touretzky, M. Mozer and M. Hasselmo (Eds). Advances in Neural\n');
fprintf('Information Processing Systems 8:145-151, MIT Press, Cambridge, MA (1996).\n');
fprintf('\n');

% Determined from the data
nsamps=size(datain.data,1);
nchans=size(datain.data,2);
ntrials=size(datain.data,3);
tlength=nsamps/datain.sr;

disp('Assuming all non-zero and non-NaN data are good...');
thevars=squeeze(var(datain.data));
artifact=thevars==0 | isnan(thevars);
goodtrials=setdiff(1:ntrials,find(sum(artifact)==nchans));
goodchans=setdiff(1:nchans,find(sum(artifact,2)==ntrials));
goodchans=setdiff(goodchans,badchans);
ngoodchans=length(goodchans);
ngoodtrials=length(goodtrials);

% Choose reasonable number of components to solve for if not given
if isempty(ncomps);
    ncomps=min([60 ngoodchans round((ngoodtrials*nsamps)^(1/3))]);
    disp(['Solving for ' num2str(ncomps) ' components...']);
end

% Keeps top 60 components. Smaller ones contribute insignificant variance
% for the most part
if isempty(nkeep);
    nkeep=min([60 ncomps]);
end

% Checks for bad inputs
if ncomps > ngoodchans;
    disp('ncomps cannot be greater than ngoodchans, resetting ncomps to ngoodchans...');
    ncomps=ngoodchans;
end
if ncomps > round((ngoodtrials*nsamps)^(1/3));
    ncomps=round((ngoodtrials*nsamps)^(1/3));
    disp(['Not enough samples per component, lowering ncomps to ' num2str(ncomps) '...']);
end
if nkeep > ncomps;
    disp('nkeep cannot be greater than ncomps, resetting nkeep to ncomps...');
    nkeep=ncomps;
end
if fftfreq>datain.sr/2;
    disp('fftfreq cannot be higher than nyquist, resetting to nyquist...');
    fftfreq=datain.sr/2;
end

% Other calculated values
nbins=ceil(fftfreq*tlength);

% Concatenate data to prepare for ICA
alldata=segtocat(datain.data(:,goodchans,goodtrials));
alldatasquared=sum(sum(alldata.^2));

% Prepare the ICA data structure
datain.ica=zeros(nsamps,ncomps,ntrials);
datain.sep=zeros(ncomps,nchans);
datain.mix=zeros(ncomps,nchans);
datain.cpvars=zeros(1,ncomps);

if extica; tmpstr='extended '; else tmpstr=''; end;
disp(['Running ' tmpstr 'Infomax ICA on the data...']);
[w,s]=runica(alldata','verbose','off','pca',ncomps,'extended',extica);

% Put ICA data into the output structure
datain.sep(:,goodchans)=w*s;
datain.mix=pinv(datain.sep)';
icasig=alldata*datain.sep(:,goodchans)';
datain.ica(:,:,goodtrials)=cattoseg(icasig,nsamps);

disp('Getting percent of variance that each component accounts for...');
% Uses the average of two different methods, one that overestimates and one
% that underestimates.
cpvars1=zeros(1,ncomps);
cpvars2=zeros(1,ncomps);
for c=1:ncomps;
    compproj=icasig(:,c)*datain.mix(c,goodchans);
    
    datalesscomp=alldata-compproj;
    sumdifsquared=sum(sum(datalesscomp.^2));
    cpvars1(c)=100*(1-sumdifsquared/alldatasquared);
    
    projsquared=sum(sum(compproj.^2));
    cpvars2(c)=100*(projsquared/alldatasquared);
end
datain.cpvars=(cpvars1+cpvars2)/2;
    
disp('Getting variance-weighted amplitude specta for each component...');
tmp=fft(datain.ica)/nsamps;
fcoefs=tmp((1:nbins)+1,:,:);

trialvars=squeeze(var(datain.ica));

tmp=zeros(nbins,ncomps,ntrials);
for c=1:ncomps;
    tmp(:,c,goodtrials)=squeeze(abs(fcoefs(:,c,goodtrials))).*(ones(nbins,1)*trialvars(c,goodtrials));
end
datain.wspecs=mean(tmp(:,:,goodtrials),3);

% Reorder based on calculated percent variance
[~,sortorder]=sort(datain.cpvars,'descend');
datain.ica=datain.ica(:,sortorder,:);
datain.sep=datain.sep(sortorder,:);
datain.mix=datain.mix(sortorder,:);
datain.cpvars=datain.cpvars(sortorder);
datain.wspecs=datain.wspecs(:,sortorder);

% Throws out components beyond nkeep
datain.ica=datain.ica(:,1:nkeep,:);
datain.sep=datain.sep(1:nkeep,:);
datain.mix=datain.mix(1:nkeep,:);
datain.cpvars=datain.cpvars(1:nkeep);
datain.wspecs=datain.wspecs(:,1:nkeep);

% Tosses orginal data 
datain=rmfield(datain,'data');