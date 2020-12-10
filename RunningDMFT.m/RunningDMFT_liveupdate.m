%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HOW TO USE?
%
%  - Put this script, together with an input-file for your driver in a path
%    > This path will contain directories for all the U values you set-up.
%  - Set-up the name of your driver program (without .f90 extension)
%    > e.g. driver = 'ed_kane_mele';
%  - Adjust Umin, Uman and Ustep to your desire --> U=Umin:Ustep:Umax
%  - Set SOI to your desire: --> you will get a fixed-SOI linear span
%  - Adjust Uold to catch a 'restart-folder' in the path [!applies -> -1]
%  - Adjust notConvThreshold to control how much DMFT will wait for you (!)
%  - Select doMPI (true.or.false) to run with openMPI or not
%  - Run everything with $ matlab -batch RunningDMFT_liveupdate
%  - Once in a while control for the presence of ERROR.README in the path
%    > If you see it DMFT is not converging; you might want to inspect the
%      the output files and update wmixing in the **external** input-file
%      to a suitable value. At the next run DMFT will catch the update ;)
%  - At the end you will find some additional output in the U=%f folders
%    > a LOG_dmft.txt which is just a mirror of the DMFT output (via tee)
%    > a LOG_time.txt which is a wall-clock-time value for the whole DMFT
%  - Also an additional output file in the main (external) path
%    > a U_list.txt that stores all the used U-values (for easier post..)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

driver = 'ed_kane_mele';	doMPI = false;

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
Umin = 1; Umax = 10;		% Input Hubbard
Ustep = 0.25;			% Hubbard updates

notConverged = 0;		% Convergence-fail *counter*
notConvThreshold = 5;		% Maximum #{times} we accept DMFT to fail
% ----------------------------> Available time to change wmix on the flight

U = Umin; Uold = -1;
doUpdate = true;
while U <= Umax                % Hubbard loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~>

UDIR= sprintf('U=%f',U);       % Make a folder named 'U=...', where '...'
mkdir(UDIR);                   % is the given value for Hubbard interaction
cd(UDIR);                      % Enter the U-folder

oldDIR=sprintf('../U=%f',Uold);% -----------------------------------------
if isfolder(oldDIR)            % If it exist a "previous" folder: 
copyfile(oldDIR);              % Copy everything from last dmft evaluation
end                            % -----------------------------------------

copyfile ../inputKANEMELE.conf % Copy inside the **external** input file

%% Run FORTRAN code (already compiled and added to PATH!) %%%%%%%%%%%%%%%%%
if doMPI
mpi = 'mpirun ';				% Control of MPI
else						% boolean flag...
mpi = [];
end
HUBBARD =sprintf(' uloc=%f',U);		% OVERRIDE of
T2 =sprintf(' t2=%f',SOI);			% PARAMETERS
OutLOG = ' | tee LOG_dmft.txt';
dmft_ed_call = [mpi,driver,HUBBARD,T2,OutLOG];
tic
system(dmft_ed_call);				% Fortran-call
chrono = toc;
file_id = fopen('LOG_time.txt','w');
fprintf(file_id,'%f\n', chrono);		% Write on time-log
fclose(file_id);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HERE WE CATCH A FAILED (unconverged) DMFT LOOP
if isfile('ERROR.README') 
    notConverged = notConverged + 1;
    delete ERROR.README
else
    fprintf(Ulist,'%f\n', U);	               % Write on U-log
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd ..                          % Exit the U-folder


if notConverged > notConvThreshold
   error('DMFT not converged: phase-span stops now!')
elseif notConverged > 0 && notConverged < notConvThreshold
   % You manage *the mixing* manually, by modifying on-the-flight
   % the **external** input-file at runtime. Very effective when feasible.
   doUpdate = false; % -------------------> So you don't want to update...
end

if doUpdate
   Uold = U;
   U = U + Ustep;              % Hubbard update  
else
   U = Uold + Ustep;		% old-Hubbard update (if nonconverged!)
end

end                            % <~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fclose(Ulist);

