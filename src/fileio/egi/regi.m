function [data rate info] = regi(file)
%
% [data rate info] = regi(file)
%
% REGI reads data from EGI (egi, hnl) .RAW files

% Copyright (C) 2006 William Winter
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

%% Record of Revisions
%   Date           Programmers               Description of change
%   ====        =================            =====================
%   2006           Bill Winter                  Original code
%   3/31/15        Michael Nunez    Fixed ability to handle multiple events
%   11/09/15       Michael Nunez       added event = [] if no event present

%Currently untest for .raw files with no events
%Please track your changes!

%%

fid = fopen(file,'rb','b');                 % open file, read header
[rate version n nevent N scale eventlabs header] = regihdr(fid); %#ok<ASGLU>
switch version                              % data type
    case 2, data = 'int16';
    case 4, data = 'float32';
end
data = fread(fid,[n + nevent,N],data);      % read file
fclose(fid);                                % close file

if ~isempty(eventlabs)
    event = struct;
    for m = 1:nevent
        event.(eventlabs(m,:)) = data(n+m,:);
    end
    %event = struct(event,data(n+1:n+nevent,:)); % separate event data
    data(n+1:n+nevent,:) = [];
else
    event = [];
end

if n == 129, data(n,:) = []; end;           % Remove 129 from NS 3.0.1 data
data(end + 1,:) = 0;                        % add reference channel
global eegchan                              % data channels
[gain zero] = regigz(file);                 % gain and zero

% data(eegchan,:)=applygz(data(eegchan,:),gain(eegchan)/400,zero(eegchan));
for k = eegchan, data(k,:) = 400*(data(k,:) - zero(k))./gain(k); end
% data = navref(data,1);                    % average reference data

info = struct('filepath',file,...           % extra information
   'version',version,'event',event,...
   'header',header);

% siz = size(data);
% out = struct('Data',data,'SamplingRate',rate,...
%     'NSample',siz(2),'NTotalSample',siz(2),...
%     'NChannel',siz(1),'ChannelLabel',cell(1,siz(1)),...
%     'ChannelScaleFactor',400./gain);