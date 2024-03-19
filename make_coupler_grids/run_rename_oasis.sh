#!/usr/bin/env bash

for x in rmp*nc ; do
    ./rename_for_oasis.py $x
done