%% Main
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

global ignUlist
ignUlist = false;

% We don't have a SOI-values list, but we can obtain that by just
% inspecting the subdirectories...
[SOI_list, SOI_names] = get_list('SOI');
Nlines = length(SOI_list);
for iSOI = 1:Nlines
    SOIDIR = SOI_names(iSOI);
    fprintf(strcat(SOIDIR));
    cd(SOIDIR);
    if isfile('U_list.txt')
        U_list = load('U_list.txt');
    else
        U_list = [];
    end
    [ids,ordpms,U_list] = extract_line(U_list); fprintf('..DONE\n');
    save('order_parameter_line.mat','ids','ordpms','U_list');
    cd('..');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [flist, strlist] = get_list(VARNAME)
%% Getting a list of variable values, from directories.
%  VARNAME: a string, identifying the listed variable (e.g. 'U')
%  flist: an array of float_values (e.g. U=[:] )
%  strlist: an array of dir_name strings (e.g. ['U=%f'] )
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subthings = dir('.'); % Retrieves every subdir and file inside pwd
    subfolders = subthings([subthings(:).isdir]); % Keeps only subfolders
    subfolders = subfolders(~ismember({subfolders(:).name},{'.','..'}));
    N = length(subfolders); flist = zeros(N,1); strlist = strings(N,1);
    for i = 1:N
        DIR = subfolders(i).name; % Let's get the indexed string...
        flist(i) = sscanf(DIR, [VARNAME,'=%f']); %...and extract the value!
        strlist(i) = DIR;
    end
    % We need to sort the lists by floats (not strings, as it is now)
    [flist, sortedIDX] = sort(flist); strlist = strlist(sortedIDX);
end

function [ids,ordpms,U_list] = extract_line(U_LIST)
%% Getting a list of variable values, from directories.
%  U_LIST: an array of values for Hubbard interaction U (could be empty!)
%  ids: a cell of strings, the names of the order parameters 
%  ordpms: a cell of arrays, corresponding to the names above, for all U
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    global ignUlist
    if isempty(U_LIST) | ignUlist == true
       [U_LIST, ~] = get_list('U'); 
    else
       U_LIST = sort(U_LIST);
    end
    % Then we can proceed spanning all the U-values
    Nu = length(U_LIST);
    cellordpms = cell(Nu,1);
    for iU = 1:length(U_LIST)
        U = U_LIST(iU);
        UDIR= sprintf('U=%f',U);
        if ~isfolder(UDIR)
           errstr = 'U_list file appears to be inconsistent: ';
           errstr = [errstr,UDIR];
           errstr = [errstr,' folder has not been found.'];
           error(errstr);
        end
        cd(UDIR); 
        [ids, cellordpms{iU}] = get_order_parameters();
        cd('..');
    end
    % We need some proper reshaping
    Nordpms = length(ids);
    ordpms = cell(1,Nordpms);
    for jORDPMS = 1:Nordpms
        ordpms{jORDPMS} = zeros(Nu,1);
        for iU = 1:Nu
           ordpms{jORDPMS}(iU) = cellordpms{iU}(jORDPMS);
        end
    end
    U_list = U_LIST;
end

function [names, order_parameters] = get_order_parameters()
%% Getting all information from order_parameters_[].dat
%  names: a cell of strings, the names of the order parameters (from [])
%  order_parameters: an array of floats, corresponding to the names above
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %names = readcell('order_parameters_info.ed','FileType','fixedwidth');
    %names(strcmp(names,'#'))=[];
    pattern = 'order_parameters_';
    files = dir('.');
    N = length(files);
    for i = 1:N
        temp_name = getfield(files,{i},'name');
        found = strfind(temp_name,pattern);
        if found
           full_name = temp_name;
        end
    end
    beheaded = erase(full_name,pattern); % Removes 'order_parameter_'
    detailed = erase(beheaded,'.dat');   % Removes '.dat'
    names = strsplit(detailed,'_');      % Reads the names separated by '_'
    order_parameters = load(full_name);
    if length(names) ~= length(order_parameters)
       error('Something went wrong reading order parameters!'); 
    end
end
