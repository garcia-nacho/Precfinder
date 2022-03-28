#!/bin/bash
if [[ "$(docker images -q precfinder 2> /dev/null)" != "" ]] 
then
  docker run -it --rm -v $(pwd)/Inference:/Inference precfinder Rscript /home/docker/Scripts/Inference.R ${1}
else
  echo "No precfinder found! Run Install.sh to install it."
fi



