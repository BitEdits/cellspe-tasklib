#!/bin/bash

aclocal
libtoolize --copy --force
autoconf
autoheader
automake -a

./configure --prefix=/usr
make clean
make
sudo make install
