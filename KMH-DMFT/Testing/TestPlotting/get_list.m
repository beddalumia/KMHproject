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