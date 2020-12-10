% Let MATLAB see the goddamn PATH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -> works only if matlab has been started from a unix terminal! (0^0~~,)
path = getenv('PATH');
path = [path ':/usr/local/bin'];
setenv('PATH', path) 
clear path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ulist = fopen('U_list.txt','w');
%% Phase-Line: single loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SOI  = 0;                      % Input Spin-Orbit 
Umin = 1; Umax = 10;           % Input Hubbard 
Ustep = 0.25;			% Phase-line step
Uold = -1;			% Restart option

U = Umin; 
while U <= Umax                % Hubbard loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~>

UDIR= sprintf('U=%f',U);       % Make a folder named 'U=...', where '...'
mkdir(UDIR);                   % is the given value for Hubbard interaction
cd(UDIR);                      % Enter the U-folder

oldDIR=sprintf('../U=%f',Uold);% -----------------------------------------
if isfolder(oldDIR)            % If it exist a "previous" folder: 
copyfile(oldDIR);              % Copy everything from last dmft evaluation
end                            % -----------------------------------------

copyfile ../inputKANEMELE.conf % Copy inside the *external* input file

%% Run FORTRAN code (already compiled and added to PATH!) %%%%%%%%%%%%%%%%%
mpi = 'mpirun ';
driver = 'ed_kane_mele';
HUBBARD =sprintf(' uloc=%f',U);		% OVERRIDE of
T2 =sprintf(' t2=%f',SOI);			% PARAMETERS
OutLOG = ' > LOG.out';
dmft_ed_call = [mpi,driver,HUBBARD,T2,OutLOG];
tic
system(dmft_ed_call);				% Fortran-call
chrono = toc;
file_id = fopen('LOG_time.txt','w');		% Write on time-log
fprintf(file_id,'%f\n', chrono);
fclose(file_id);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HERE WE NEED TO CATCH SOMEHOW A FAILED (unconverged) DMFT LOOP
if isfile('ERROR.README') 
   fclose(Ulist);
   error('DMFT not converged: phase-span stops now!')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd ..                          % Exit the U-folder

fprintf(Ulist,'%f\n', U);	% Write on U-log

Uold = U;
U = U + Ustep;              	% Hubbard update  

end                            % <~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fclose(Ulist);

