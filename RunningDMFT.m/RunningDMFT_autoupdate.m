% Let MATLAB see the goddamn PATH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -> works only if matlab has been started from a unix terminal! (0^0~~,)
path = getenv('PATH');
path = [path ':/usr/local/bin'];
setenv('PATH', path) 
clear path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Phase Line: single loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

copyfile ../inputKANEMELE.conf % Copy inside the *external* input file

%% Run FORTRAN code (already compiled and added to PATH!) %%%%%%%%%%%%%%%%%
mpi = 'mpirun ';
driver = 'ed_kane_mele';
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
%% HERE WE NEED TO CATCH SOMEHOW A FAILED (unconverged) DMFT LOOP
if isfile('ERROR.README')
    notConverged = notConverged + 1;
    if Uold >= Umin
    U = Uold; % if nonconverged we don't want to update! 
    end
    delete ERROR.README
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
