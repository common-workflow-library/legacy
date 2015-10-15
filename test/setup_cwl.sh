#!/usr/bin/env bash

pwd
echo $PATH
echo $*

. /vagrant/cwl.conf

cd "$BASEDIR"

setup_cwl() {
cd reference && sudo python setup.py install
cd cwl-runner && sudo python setup.py install
}

cd "${DIST_DIR}" && git pull || { git clone "${GIT_URL}" "${DIST_DIR}" && cd "${DIST_DIR}"; }
setup_cwl

cd "$BASEDIR"

