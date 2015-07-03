#!/bin/sh

# --force : rebuild configure script no matter timestamp related to configure.ac
# --install : copy some missing files (e.g. COPYING, INSTALL)
# -I config -I m4 : files are placed in those subdirectories
mkdir -p config m4
autoreconf --force --install -I config -I m4
