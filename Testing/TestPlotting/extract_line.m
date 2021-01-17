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
        [ids, obs{iU}] = get_observables();
        cd('..');
    end
    U_list = U_LIST;
end