HERE=$(pwd) # Saving in variable HERE the current directory path (pwd: print working directory)
while read t2; #<-------------------------------------------------------------------------------
do
  T2DIR=SOI=$t2;
  mkdir $T2DIR;
  cp list_u inputKANEMELE.conf $T2DIR;

  while read u; #<~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  do
    UDIR=U=$u;
    mkdir $MDIR/$UDIR;
    cp inputKANEMELE.conf *.restart $T2DIR/$UDIR/;
    cd $T2DIR/$UDIR;
    echo $t2 $u;
    echo "Executing: ed_kane_mele t2=$t2 uloc=$u |tee log.out"
    ed_kane_mele t2=$t2 uloc=$u |tee log.out
    cd $HERE;
  done <list_u #>~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

done <list_t2  #>--------------------------------------------------------------------------------

# Formatting of the lists:
# both list_u and list_t2 are just one-column vectors of real values (Hubbard U and SOI, resp.)
