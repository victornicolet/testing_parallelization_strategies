#!/usr/bin/env bash

for numcore in 1 2 3 4 5 6 7 8 9 10 12 14 16 18 20 24 28 32; do
  ./mbo_test 28 numcore
done