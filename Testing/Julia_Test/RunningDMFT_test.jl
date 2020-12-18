# run(`ed_kane_mele | tee LOG.dmft`) 
# -> DOES NOT WORK: 
#    ERROR: LoadError: LoadError: parsing command `ed_kane_mele | tee LOG.dmft`:
#    special characters "#{}()[]<>|&*?~;" must be quoted in commands

# run(`ed_kane_mele \| tee LOG.dmft`)
# -> DOES NOT WORK: 
#    Runs and all, but LOG.dmft never happens to exist...

# run(`ed_kane_mele \> tee LOG.dmft`)
# -> DOES NOT WORK: 
#    Runs and all, but LOG.dmft never happens to exist...also stdout actually is printed on screen.

# run(`ed_kane_mele > tee LOG.dmft`)
# -> DOES NOT WORK: 
#    ERROR: LoadError: LoadError: parsing command `ed_kane_mele | tee LOG.dmft`:
#    special characters "#{}()[]<>|&*?~;" must be quoted in commands

# run(`bash ed_kane_mele \| tee LOG.dmft`)
# -> DOES NOT WORK: 
#    ERROR: LoadError: failed process: 
#    Process(`bash ed_kane_mele '|' tee LOG.dmft`, ProcessExited(126)) [126]

run(`bash -c "ed_kane_mele | tee LOG.dmft" `)
# -> Finally, this works.

