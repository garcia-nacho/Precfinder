#!/bin/bash
if [[ "$(docker images -q garcianacho/precfinder 2> /dev/null)" != "" ]] 
then
  docker run -it --rm -v Inference:/Inference garcianacho/recfinder Rscript /home/docker/Scripts/Inference.R ${1}
else
  echo "No precfinder found! Run Install.sh to install it."
fi



