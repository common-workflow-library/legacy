#!/usr/bin/env bash

. /vagrant/cwl.conf

cd "${BASEDIR}/${DIST_DIR}/${TEST_DIR}"

./run_test.sh CWLTOOL=../reference 
