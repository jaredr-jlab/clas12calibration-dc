#!/bin/sh -f

OS=$(uname)

# set up environment
classPath="$PWD/target/*:$PWD/target/classes/*"

# run 
java -Xmx1536m -Xms1024m  -Dsun.java2d.pmoffscreen=false -Djava.util.logging.config.file=$CLAS12DIR/etc/logging/debug.properties -cp $classPath org.clas.detector.clas12calibration.viewer.Driver
