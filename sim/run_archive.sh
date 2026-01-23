#!/bin/bash

REPO_NAME=halld_sim_srcct
GIT_TAG=5.1.0.2
BUILD_TAG=hdr512

echo cd /work/halld2/home/boyu/src_software_builds/${REPO_NAME}
if [[ $BUILD_TAG == "none" ]]; then
    echo git archive --prefix=${REPO_NAME}-${GIT_TAG}/ -o ../temp.tar.gz ${GIT_TAG}_srcct
else
    echo git archive --prefix=${REPO_NAME}-${GIT_TAG}^${BUILD_TAG}/ -o ../temp.tar.gz ${GIT_TAG}_srcct
fi

echo cd ..
echo tar -xvzf temp.tar.gz
echo rm temp.tar.gz