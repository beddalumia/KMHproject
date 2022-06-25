%% HOW TO USE?
%
%  - Put this script, together with an input-file for your model, in a path
%    > This path will contain directories for all the U values you set-up.
%  - Set-up the name of your driver program (without .f90 extension)
%    > e.g. driver = 'mf_km_2d';
%  - Adjust Umin, Umax and Ustep to your desire --> U = Umin:Ustep:Umax
%  - Set SOI to your desire: --> you will get a fixed-SOI linear span
%  - Adjust Uold to catch a 'restart-folder' in the path [!applies -> -1]
%  - Select doMPI (true.or.false) to run with openMPI or not
%  - Run everything with $ matlab -batch DryMF_line
%  - At the end you will find some additional output in the U=%f folders
%    > a LOG_out.txt which is just a mirror of the MF output (via tee)
%    > a LOG_time.txt which is a wall-clock-time value for the whole MF
%  - Also an additional output file in the main (external) path
%    > a U_list.txt that stores all the used U-values (for easier post..)
%      NB. Only the converged calculations will update the U_list :D
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

driver  = 'mf_km_afm';
doMPI   = false;        % >> MF-code is not MPI-safe <<

SOI     = 0.00;         % Fixed Spin-Orbit 
Uold    = 2.00;			% Restart point
Umin    = 2.05;         % Phase-line start
Ustep   = 0.05;			% Phase-line step
Umax    = 4.00;         % Phase-line end
Nopms   = 4;            % #{order parameters}

runDMFT.dry_line(driver,doMPI,Uold,Umin,Ustep,Umax,'t2',SOI,'nparams',Nopms);
