#!/bin/bash

#conda update pangolin
pangolin --update
pangolin ${1} -t 10 --outfile ${1}_pango.csv


