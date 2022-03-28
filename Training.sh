#!/bin/bash

if [[ "$(docker images -q garcianacho/precfinder 2> /dev/null)" != "" ]] 
then
  mkdir {1}
  docker run -it --rm -v $(pwd):/Inference -v $(pwd):/Models -v $(pwd)/Training/:/Training garcianacho/recfinder Rscript /home/docker/Scripts/Training.R {1}
else
    echo "Not docker image found, build docker image first!"
fi
