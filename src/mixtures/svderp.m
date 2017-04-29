%%%%%%%%%%UNFINISHED%%%%%%%%%%%
% function newdata = svderp(datain,varargin)
%
% This function performs an SVD decomposition of segmented data.
%
% The returned data structure has all of the original fields, plus now
% includes unmixing matrices, as well as amplitude spectra
% for each component and the percentage of variance that each component
% accounts for.
%
% Required Input:
%    Data structure containing the fields .data and .sr 
%
% Optional Arguments: 
%
%    erpwind - the window to average the data (calculate the ERP).
%              Default: [1:size(datain.data,1)]
%
%    baseline - window to calculate the baseline of the ERP
%               Default: [] (no baseline)
%
%    ncomps-  the number of components to solve for. Default: nchans
%
%    fftfreq- the max frequency to calculate for the component amplitude
%             spectra.  Default: 50
%
%    badchans - channels that are not zeroed out but should be considered
%               bad and excluded from the ICA.
%                    
%
% Copyright (C) 2017 Michael D. Nunez, <mdnunez1@uci.edu>
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

function datain = svderp(datain,varargin)

%% Record of Revisions
%   Date           Programmers               Description of change
%   ====        =================            =====================
%  04/28/17       Michael Nunez          Converted from icasegdata.m


if nargin < 1; help svderp; return; end;

% Parse inputs;
[~,erpwind,baseline,ncomps,fftfreq,badchans]=...
    parsevar(varargin,'erpwind',1:size(datain.data,1),...
        'baseline',[],'ncomps',[],...
        'fftfreq',50,'badchans',[]);

fprintf('SVD ERP used! Please cite:\n');
fprintf('\n');
fprintf('Nunez, M. D., Vandekerckhove, J., & Srinivasan, R., \n');
fprintf('How attention influences perceptual decision making:\n');
fprintf('Single-trial EEG correlates of drift-diffusion model parameters.\n');
fprintf('Journal of Mathematical Psychology, 76, 117-130. (2017).\n');
fprintf('\n');


% Determined from the data
nsamps=size(datain.data,1);
nchans=size(datain.data,2);
ntrials=size(datain.data,3);
tlength=nsamps/datain.sr;

disp('Assuming all non-zero and non-NaN data are good...');
thevars=squeeze(var(datain.data));
artifact=thevars==0 | isnan(thevars);

if isfield(datain,'artifact')
     disp('Rejecting data from .artifact field matrix...');
     artifact = datain.artifact | artifact;
end

goodtrials=setdiff(1:ntrials,find(sum(artifact)==nchans));
goodchans=setdiff(1:nchans,find(sum(artifact,2)==ntrials));
goodchans=setdiff(goodchans,badchans);
ngoodchans=length(goodchans);
ngoodtrials=length(goodtrials);

% Choose reasonable number of components to solve for if not given
if isempty(ncomps);
    ncomps=ngoodchans;
    disp(['Solving for ' num2str(ncomps) ' components...']);
end

% Smaller components contribute insignificant variance for the most part
if isempty(nkeep);
    nkeep=ncomps;
end

% Checks for bad inputs
if ncomps > ngoodchans;
    fprintf('ncomps cannot be greater than ngoodchans, resetting ncomps to %d...\n',ngoodchans);
    ncomps=ngoodchans;
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