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
%    erptrials - Trial index to calculate erp
%                Default: All non-artifact trials
%
%    ncomps-  the number of components to solve for. Default: nchans
%
%    fftfreq- the max frequency to calculate for the component amplitude
%             spectra.  Default: 50
%
%    badchans - channels that are not zeroed out but should be considered
%               bad and excluded from the SVDERP.
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
%  05/01/17       Michael Nunez          Converted from icasegdata.m


if nargin < 1; help svderp; return; end;

% Parse inputs;
[~,erpwind,baseline,erptrials,ncomps,fftfreq,badchans]=...
    parsevar(varargin,'erpwind',1:size(datain.data,1),...
        'baseline',[],'erptrials',1:size(datain.data,3),'ncomps',[],...
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

% Prepare the SVDERP data structure
datain.cmp=zeros(nsamps,ncomps,ntrials);
datain.sep=zeros(ncomps,nchans);
datain.cpvars=zeros(1,ncomps);

%Calculate the ERP
erptrials = intersect(goodtrials,erptrials);
dataforerp = datain.data(erpwind,:,erptrials);
dataforerp(:, setdiff(1:nchans,goodchans), :) = 0;
if ~isempty(baseline)
    dataforerp = eegBaseline(dataforerp,baseline);
end
erp = mean(dataforerp,3);

%Use singular value decomposition
[U,S,V] = svd(erp,'econ');

% Concatenate data to prepare for multiplication
alldata=segtocat(datain.data);

% Put SVDERP data into the output structure
datain.sep=V(:, 1:ncomps)';
datain.mix=datain.sep;
cmpsig=alldata*datain.sep';
datain.cmp=cattoseg(cmpsig,nsamps);

disp('Getting percent of variance that each component accounts for...');
cpvars=diag(S.^2)/sum(diag(S.^2))*100;
datain.cpvars = cpvars(1:ncomps);
    
disp('Getting variance-weighted amplitude specta for each component...');
tmp=fft(datain.cmp)/nsamps;
fcoefs=tmp((1:nbins)+1,:,:);

trialvars=squeeze(var(datain.cmp));

tmp=zeros(nbins,ncomps,ntrials);
for c=1:ncomps;
    tmp(:,c,goodtrials)=squeeze(abs(fcoefs(:,c,goodtrials))).*(ones(nbins,1)*trialvars(c,goodtrials));
end
datain.wspecs=mean(tmp(:,:,goodtrials),3);

% Reorder based on calculated percent variance
[~,sortorder]=sort(datain.cpvars,'descend');
datain.cmp=datain.cmp(:,sortorder,:);
datain.sep=datain.sep(sortorder,:);
datain.mix=datain.mix(sortorder,:);
datain.cpvars=datain.cpvars(sortorder);
datain.wspecs=datain.wspecs(:,sortorder);

% Tosses orginal data 
datain=rmfield(datain,'data');