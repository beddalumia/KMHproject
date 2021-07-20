clear all
clc
U = 7; UID = sprintf('U=%f',U);
observables = load('observables_all.ed');
Nloops = max(size(observables));
loops = 1:Nloops;
        %%% Structure (for ed_kane_mele) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        name{1} =  'Density';                
        name{2} =  'Double Occupancy'; 
        name{3} =  'n_{up}';                           
        name{4} =  'n_{dw}';                              
        name{5} =  'n_{up}-n_{dw}';                                
        name{6} =  '<S^2>';
        name{7} =  'E_{gs}';
        name{8} =  '<Sz^2>';
        name{9} =  'N2_11';  
        name{10}=  'Z_{up}'; 
        name{11}=  'Z_{dw}';           
        name{12}=  'ImSigma_{up}(0)';
        name{13}=  'ImSigma_{dw}(0)';
        name{14}=  'Nph';
        name{15}=  'Wph';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for id = [1,2,5,8,10]
    fig = figure("Name",name{id},'visible','off');
    plot(loops,observables(:,id));
    xlabel('loops');
    ylabel(name{id});
    title(UID);
    fprintf('Printing %s\n', name{id});
    print(append(UID,'_',name{id},'.png'),'-dpng','-r600')
    close(fig);  
end
