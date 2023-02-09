#!/bin/bash
if [ -z "$1" ];
then
  echo "No model name selected. Run it as Training.sh YourModelName" 
else

if [[ "$(docker images -q precfinder 2> /dev/null)" != "" ]] 
then
  mkdir Models/${1}
  docker run -it --rm -v $(pwd)/Models:/Models -v $(pwd)/Training:/Training precfinder Rscript /home/docker/Scripts/Training.R ${1}
  echo "Adding your model to the image..."
  docker build -t precfinder .
else
   echo "No precfinder found! Run Install.sh to install it."
fi

fi
