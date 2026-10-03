# Pseudopotentials

The PRA manuscript specifies Troullier-Martins norm-conserving pseudopotentials
in Kleinman-Bylander form. Set the actual Li and F filenames under
`pseudopotentials` in the campaign configuration and supply those files in
each run's `input` directory. No replacement files or unverified filenames
are supplied. The runner records their checksums; the declared family must
also be checked against the actual source files before a physical calculation.
