#!/usr/bin/bash

LD_LIBRARY_PATH=${HOME}/.xacc/lib:${HOME}/evendim/src/xacc  ./testXacc --accelerator-name tnqvm
