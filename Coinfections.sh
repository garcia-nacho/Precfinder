#!/bin/bash
if [[ "$(docker images -q precfinder 2> /dev/null)" != "" ]] 
then

if [[ ${1} == "Nanopore" ]]
then
  docker run -it --rm -v $(pwd)/Coinfections:/Noise precfinder Rscript /home/docker/Scripts/NoiseNanopore.R
fi

if [[ ${1} == "Illumina" ]]
then
  docker run -it --rm -v $(pwd)/Coinfections:/Noise precfinder Rscript /home/docker/Scripts/NoiseIllumina.R
fi

if [ ${1} != "Nanopore" ]  && [ ${1} != "Illumina" ]
then
  echo "Try again with Nanopore or Illumina as argument (i.e. Coinfections.sh Illumina or Coinfections.sh Nanopore)"
fi

else
  echo "No precfinder found! Run Install.sh to install it."
fi
