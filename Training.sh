#!/bin/bash

docker run -it --rm -v $(pwd)/Inference:/Inference -v $(pwd)/Models/:/Models -v $(pwd)/Training/:/Training garcianacho/recfinder Rscript /home/docker/Scripts/Training.R
