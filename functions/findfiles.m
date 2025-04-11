
% findfiles: Locates BLE data files and maps them to station indices
%   Searches specified folder for JSON files matching BLE packet naming
%   convention, extracts station IDs, and maps them to predefined station list
%
%
%   Coded by Yi Jin

%   Inputs:
%   - folderPath [string]
%       Path to directory containing BLE IQ data files
%   - id_Stations [cell/string array]
%       List of known station identifiers
%   - pstr [string]
%       Point string for file filtering
%
%   Outputs:
%   - matchingFiles_s [cell array]
%       Sorted list of matching filenames
%   - stationIndices_s [array]
%       Corresponding station indices in original id_Stations order

function [ matchingFiles_s, stationIndices_s ] = findfiles( folderPath, id_Stations, pstr  )

% Step 1: Find files
allFiles = dir(folderPath);
matchingFiles = {};


for k = 1:length(allFiles)
    fileName = allFiles(k).name;
    
    if contains(fileName, ['IQ_',pstr,'_ble'])
        matchingFiles{end+1} = fileName;
    end
end

% Step 2: Extract station numbers
stationIndices = zeros(1, length(matchingFiles));
matchingFiles_s = cell( length(matchingFiles),1 ); 
%
for i = 1:length(matchingFiles)
    filename = matchingFiles{i};
    id_start = strfind(filename, '_ble-pd-') + 1;
    id_end = strfind(filename, '.json') - 1;
    stationID = filename(id_start:id_end);
    
    % Find ID
    index = find(contains(cellstr(id_Stations), stationID));
    
    % save results
    if ~isempty(index)
        stationIndices(i) = index;
    else
        stationIndices(i) = NaN; 
    end
end


[stationIndices_s, idx] = sort(stationIndices);
for ii = 1:length(stationIndices)
    matchingFiles_s{ii} = matchingFiles{idx(ii)};
end


end



