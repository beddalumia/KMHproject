id = fopen('LOG.time','wt');
tic
%% Define on-the-flight updates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = 0;
SOI = 0;
Wmix = 0.1;
%% Run FORTRAN code (already compiled and added to PATH!) %%%%%%%%%%%%%%%%%
driver = 'ed_kane_mele';
HUBBARD =sprintf(' uloc=%f',U);		% OVERRIDE
T2 =sprintf(' t2=%f',SOI);			% of
MIXING = sprintf(' wmixing=%f',Wmix);		% PARAMETERS
OutLOG = ' > LOG.dmft';
dmft_ed_call = [driver,HUBBARD,T2,MIXING,OutLOG];
ed_status = system(dmft_ed_call)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ETA = toc;
fprintf(id,'%f\n', ETA);
fclose(id);
