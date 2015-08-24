#!/bin/bash

export LD_LIBRARY_PATH=../../:$LD_LIBRARY_PATH
./bounds data/n6_all.txt swapBip 1e-3
