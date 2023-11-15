#!/usr/bin/bash

LD_LIBRARY_PATH=${HOME}/.xacc/lib:${HOME}/evendim/src/xacc  ./QuantumGEPXacc  --accelerator-name tnqvm
