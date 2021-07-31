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