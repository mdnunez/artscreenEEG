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
%             data. Default: min([nchans round(goodsamps^(1/3))])
%
%    nkeep-   the number of components to retain for review. The comps
%             are sorted by variance, thus only minor artifacts comprise
%             the later components. Defaut: ncomps
%
%    fftfreq- the max frequency to calculate for the component amplitude
%             spectra.  Default: 50
%
%    algorithm - string indicating which ICA algorithm to use. 
%                Choices: 'infomax' or 'fastica'. Default: 'fastica'
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
%   1.4 - Changed the default number of components to solve for. 
%         8/23/16 -Michael Nunez 
%   1.5 - Export citation to the terminal when Infomax ICA. 
%         12/20/16 - Michael Nunez
%   1.6 - Change defaults for the number of components to solve for 
%         and keep, reject from datain.artifact 04/24/17 - Michael Nunez
%   1.7 - Addition of FastICA algorithm 05/02/17 - Michael Nunez
%   1.8 - Change default to FastICA algorithm, adding try statement for 
%         Infomax ICA 5/15/17 - Michael Nunez
%   1.9 - Fixing FastICA ncomp and nkeep error 7/11/17 - Michael Nunez
%

%To do:
% 1) Fix Infomax ICA matrix multiplication
% 2) Verify use of Moore-Penrose pseudoinverse after Infomax ICA
% 3) Add SOBI as a feature to artscreenEEG

if nargin < 1; help icasegdata; return; end;

% Parse inputs;
[~,ncomps,nkeep,fftfreq,algorithm,extica,badchans,verbose]=...
    parsevar(varargin,'ncomps',[],'nkeep',[],...
    'fftfreq',50,'algorithm','fastica','extica',1,'badchans',[],'verbose','off');

if (~strcmp(algorithm,'infomax')) & (~strcmp(algorithm,'fastica'))
    error('Algorithm choices are either ''infomax'' or ''fastica''');
end

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
    ncomps=min([ngoodchans round((ngoodtrials*nsamps)^(1/3))]);
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
if ncomps > round((ngoodtrials*nsamps)^(1/3));
    suggested_ncomps=round((ngoodtrials*nsamps)^(1/3));
    warning(['Possibly not enough samples per component. It is recommended to lower ncomps to ' num2str(suggested_ncomps) '...']);
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

infomaxerr = 0;
if strcmp(algorithm,'infomax')
    try
        if extica; tmpstr='extended '; else tmpstr=''; end;
        disp(['Running ' tmpstr 'Infomax ICA on the data...']);
        [w,s]=runica(alldata','verbose',verbose,'pca',ncomps,'extended',extica);
        % Put ICA data into the output structure
        datain.sep(:,goodchans)=w*s;
        datain.mix=pinv(datain.sep)'; %Check use of Moore-Penrose pseudoinverse
        icasig=alldata*datain.sep(:,goodchans)';
        datain.ica(:,:,goodtrials)=cattoseg(icasig,nsamps);

        fprintf('Infomax ICA used! Please cite:\n');
        fprintf('\n');
        fprintf('Makeig, S., Bell, A.J., Jung, T-P and Sejnowski, T.J., \n');
        fprintf('Independent component analysis of electroencephalographic data.\n');
        fprintf('In: D. Touretzky, M. Mozer and M. Hasselmo (Eds). Advances in Neural\n');
        fprintf('Information Processing Systems 8:145-151, MIT Press, Cambridge, MA (1996).\n');
        fprintf('\n');
    catch me
        rethrow(me);
        fprintf('Error given by Infomax ICA. Switching to FastICA algorithm...\n');
        infomaxerr = 1;
    end
end
if strcmp(algorithm,'fastica') || infomaxerr == 1
    disp(['Running FastICA on the data...']);
        [outsig, A, W]=fastica(alldata','verbose',verbose,'lastEig', ncomps, 'numOfIC', ncomps);
        if ncomps ~= size(W,1);
            fprintf('FastICA only found %d components, resetting ncomps to %d...\n',ncomps,ncomps);
            ncomps = size(W,1);
            if nkeep > ncomps;
                disp('nkeep cannot be greater than ncomps, resetting nkeep to ncomps...');
                nkeep=ncomps;
            end
            datain.ica=zeros(nsamps,ncomps,ntrials);
            datain.sep=zeros(ncomps,nchans);
            datain.mix=zeros(ncomps,nchans);
            datain.cpvars=zeros(1,ncomps);
        end
        % Put ICA data into the output structure
        datain.sep(:,goodchans)=W;
        datain.mix(:,goodchans)=A';
        icasig = outsig';
        datain.ica(:,:,goodtrials)=cattoseg(icasig,nsamps);

        fprintf('FastICA used! Please cite:\n');
        fprintf('\n');
        fprintf('Oja, E., & Hyvarinen, A., \n'); 
        fprintf('A fast fixed-point algorithm for independent component analysis.\n');
        fprintf('Neural computation, 9(7), 1483-1492. (1997).\n');
        fprintf('\n');
end



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