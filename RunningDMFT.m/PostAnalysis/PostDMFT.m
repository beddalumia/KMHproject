%% Main
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

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
    [ids,obs,U_list] = extract_line(U_list); fprintf('..DONE\n');
    save('observables_line.mat','ids','obs','U_list');
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

function [ids,obs,U_list] = extract_line(U_LIST)
%% Getting a list of variable values, from directories.
%  U_LIST: an array of values for Hubbard interaction U (could be empty!)
%  ids: a cell of strings, the QcmPlab names of the observables 
%  obs: a cell of float-arrays, corresponding to the names above, forall U
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if isempty(U_LIST) % We shall make it, looking at subdirectories...
       [U_LIST, ~] = get_list('U'); 
    end
    % Then we can proceed spanning all the U-values
    Nu = length(U_LIST);
    cellobs = cell(Nu,1);
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
        [ids, cellobs{iU}] = get_observables();
        cd('..');
    end
    % We need some proper reshaping
    Nobs = length(ids);
    obs = cell(1,Nobs);
    for jOBS = 1:Nobs
        obs{jOBS} = zeros(Nu,1);
        for iU = 1:Nu
           obs{jOBS}(iU) = cellobs{iU}(jOBS);
        end
    end
    U_list = U_LIST;
end

function [names, observables] = get_observables()
%% Getting all information from observables_last.ed and observables_info.ed
%  names: a cell of strings, the QcmPlab names of the observables
%  observables: an array of float values, corresponding to the names above
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    names = readcell('observables_info.ed','FileType','fixedwidth');
    names(strcmp(names,'#'))=[];
    for i = 1:length(names)
        tempstr = names{i};                 % Temporary string variable
        head = sscanf(tempstr,'%d');        % Extracts the initial integer
        head = int2str(head);               % Int to Str conversion
        names{i} = erase(tempstr,head);     % Proper beheading ;D
    end
    observables = load('observables_last.ed');
end
