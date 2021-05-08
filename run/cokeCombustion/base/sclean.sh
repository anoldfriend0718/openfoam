#!/bin/bash
cd ${0%/*} || exit 1    # Run from this directory

foamCleanTutorials
rm postProcess -rf
rm *.log