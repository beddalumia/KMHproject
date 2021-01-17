function [names, observables] = get_observables()
%% Getting all information from observables_last.ed and observables_info.ed
%  names: a cell of strings, the QcmPlab names of the observables
%  observables: an array of float values, corresponding to the names above
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    names = readcell('observables_info.ed','FileType','fixedwidth');
    names(strcmp(names,'#'))=[];
    observables = load('observables_last.ed');
end