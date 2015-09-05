#!/usr/bin/env bash

pwd
echo $PATH
echo $*

. /vagrant/rabix.conf

cd "$BASEDIR"

setup_rubix() {
sudo python setup.py install
}

cd "${DIST_DIR}" && git pull || { git clone "${GIT_URL}" "${DIST_DIR}" && cd "${DIST_DIR}"; }
setup_rubix

cd "$BASEDIR"

