#!/bin/bash
if [[ "$(docker images -q garcianacho/precfinder 2> /dev/null)" != "" ]] 
then
  docker run -it --rm -v $(pwd):/Inference garcianacho/recfinder Rscript /home/docker/Scripts/Inference.R ${1}
else
    echo "Not docker image found, build docker image first!"
fi



