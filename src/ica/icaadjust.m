% function dataout = icaadjust(datain)
%
% This function is used to automatically review ICA components generated 
% by icasegdata.m and decide if each is artifact or brain-generated data.
% This algorithm uses the ADJUST algorithm.
%
% See:
%
% Mognon A, Bruzzone L, Jovicich J, Buiatti M.
% ADJUST: An Automatic EEG artifact Detector based on the Joint Use of Spatial and Temporal features.
% Psychophysiology 48 (2), 229-240 (2011).
%
% Automatic Part 3 of artscreenEEG's basic data cleaning functions:
%    artscreen.m => icasegdata.m => icaadjust.m => icatochan.m (optional)
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


function datain = icaadjust(datain)

% VERSION HISTORY:
%     1.0 - created by Michael Nunez - 8/14/17


if nargin < 1; help icaadjust; return; end;

if isfield(datain,'data') && ~isfield(datain,'ica');
    datain.ica=datain.data;
    datain=rmfield(datain,'data');
end

% Get required variables
ncomps=size(datain.ica,2);

% If no channel positions, return error
if ~isfield(datain,'hm');
    error('No ''hm'' field found. ADJUST algorithm requires a head model!\n');
end

% Prepare for Component Review
if ~isfield(datain,'compevals')
    datain.compevals=NaN(1,ncomps);
end



fprintf('ADJUST used! Please cite:\n');
fprintf('Mognon A, Bruzzone L, Jovicich J, Buiatti M, \n');
fprintf('ADJUST: An Automatic EEG artifact Detector based on the Joint Use of Spatial and Temporal features.\n');
fprintf('Psychophysiology 48 (2), 229-240 (2011).\n');
EEG = toADJUST(datain);
art = ADJUST(EEG);
datain.compevals(art)=0;


% Nested Functions*************************************************************
    function EEG = toADJUST(datain)
        % Converts data structure to one that can be read by the ADJUST1.1 algorithms
        EEG.data = datain.ica(1:2,1:2,:); %Note that this copy isn't really necessary as AJDUST just uses it to calculate data size
        EEG.icaact = permute(datain.ica,[2,1,3]);
        EEG.icawinv = datain.mix';
        nchans = size(datain.mix,2);
        %Note that this is a stupid way to organize electrode positions, and ADJUST just converts these back into a matrix
        %%MATLAB/EEGLAB likes EEG Cartesian coordinates in the following framework:
        %-positive X is towards the nose
        %-positive Y is towards the left ear
        %-positive Z is towards the vertex
        %%artscreen/EGI coordinates have been placed in the following framework (check this for new head models):
        %%figure; scatter3(datain.hm.Electrode.CoordOnSphere(:,1),datain.hm.Electrode.CoordOnSphere(:,2),datain.hm.Electrode.CoordOnSphere(:,3));
        %-positive X is towards the right ear
        %-positive Y is towards the nose
        %-positive Z is towards the vertex
        %The following changes reflect these differences
        for n=1:nchans,
            EEG.chanlocs(1,n) = struct('X',datain.hm.Electrode.CoordOnSphere(n,2),...
                'Y',-datain.hm.Electrode.CoordOnSphere(n,1),...
                'Z',datain.hm.Electrode.CoordOnSphere(n,3));
        end
        [temptheta, tempphi, tempradius] = cart2sph([EEG.chanlocs(1,:).X],...
            [EEG.chanlocs(1,:).Y],[EEG.chanlocs(1,:).Z]);
        for n=1:nchans,
            % EEG.chanlocs(1,n).sph_theta = (180/pi)*temptheta(n);
            % EEG.chanlocs(1,n).sph_phi = (180/pi)*tempphi(n);
            % EEG.chanlocs(1,n).sph_radius = tempradius(n);
            %The following was obtained from EEGLAB's convertlocs.m, line 152, and EEGLAB's sph2topo.m, lines 92 and 93
            EEG.chanlocs(1,n).theta = -temptheta(n)*(180/pi);
            EEG.chanlocs(1,n).radius = 0.5 - (tempphi(n)*(180/pi))/180;
        end
    end % end of toADJUST

end % end of main function
