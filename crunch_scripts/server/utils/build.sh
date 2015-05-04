#!/bin/bash

gcc -I /usr/include/python2.7 -c -fpic *.c
gcc -shared -o bitset.so *.o
