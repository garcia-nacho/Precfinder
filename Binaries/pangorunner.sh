#!/bin/bash


conda update pangolin
pangolin ${1} -t 10 --outfile ${1}_pango.csv


