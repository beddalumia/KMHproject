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
%    > If you see it, DMFT is not converging; you might want to inspect the
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
Umin = 0; Umax = 10;           % Input Hubbard 
Wmix = 0.3;		 	% Input Self-Mixing

Ustep = [0.5, 0.25, 0.1, 0.05, 0.01]; % To be auto-determined...hopefully!

notConverged = 0;		% Convergence-fail *counter*
notConvThreshold = 4;          % Maximum #{times} we accept DMFT to fail

U = Umin; Uold = -1;
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
HUBBARD =sprintf(' uloc=%f',U);		% OVERRIDE
T2 =sprintf(' t2=%f',SOI);			% of
MIXING = sprintf(' wmixing=%f',Wmix);		% PARAMETERS
outLOG = ' > LOG_dmft.txt';
dmft_ed_call = [mpi,driver,HUBBARD,T2,MIXING,outLOG];
tic
system(dmft_ed_call);
chrono = toc;
file_id = fopen('LOG_time.txt','w');
fprintf(file_id,'%f\n', chrono);
fclose(file_id);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HERE WE CATCH A FAILED (unconverged) DMFT LOOP
if isfile('ERROR.README')
    notConverged = notConverged + 1;
    if Uold >= Umin
    U = Uold; % if nonconverged we don't want to update! 
    end
    movefile('ERROR.README','../ERROR.README');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd ..                          % Exit the U-folder


if notConverged > notConvThreshold
   error('DMFT not converged: phase-span stops now!');         
else
   Uold = U; % if converged we update Uold and proceed to...         
end

U = U + Ustep(notConverged+1); % ...Hubbard update  

end                            % <~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fclose(Ulist);
