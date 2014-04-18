grep -IR --exclude-dir={external,scripts,\.git} '\bint\b' .
grep -IR --exclude-dir={external,scripts,\.git} '\bdouble\b' .
grep -IR --exclude-dir={external,scripts,\.git} '\bfloat\b' .
grep -IR --exclude-dir={external,scripts,\.git} 'MPI_INT' .
grep -IR --exclude-dir={external,scripts,\.git} 'MPI_DOUBLE' .
grep -IR --exclude-dir={external,scripts,\.git} 'MPI_FLOAT' .
