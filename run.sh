#!/bin/bash

################################################################################
#                                                                              #
# TRexFitter run                                                               #
#                                                                              #
################################################################################
#                                                                              #
# LICENCE INFORMATION                                                          #
#                                                                              #
# This program executes a simple run of TRexFitter.                            #
#                                                                              #
# 2016 Will Breaden Madden, w.bm@cern.ch                                       #
#                                                                              #
# Subject to ATLAS Data Access Policy, this software is released under the     #
# terms of the GNU General Public License version 3 (GPLv3).                   #
#                                                                              #
# For a copy of the ATLAS Data Access Policy, see                              #
# DOI: 10.7483/OPENDATA.ATLAS.T9YR.Y7MZ or http://opendata.cern.ch/record/413. #
#                                                                              #
# For a copy of the GNU General Public License, see                            #
# http://www.gnu.org/licenses/.                                                #
#                                                                              #
################################################################################

name="TRexFitter run"
version="2016-06-28T1731Z"

if [ -z "${1}" ]; then
    run_name="ttHbb"
else
    run_name="${1}"
fi

clear
cat << 'EOF'
_________________________________________________________________________________________________________
| _____________________________________________________________________________________________________  |
| |                                              .--==-````=-.                                         | |
| |                                            .`      (`o',  "===-..                                  | |
| |                                           /         '"`          ',                                | |
| |                                  _....__.`                       \ :                               | |
| |                             .-=-`              (                   |                               | |
| |                    _.-==---`               \    \   _.-.vv,        :                               | |
| |                  .`                         \    ';`    ^ VVV\/\/\/                                | |
| |               _.`              ,             \    '!    \^   ``                                    | |
|.-""""'--,,__,,-"                  /      ,      .,   '!    \^        _____________________________   | |
/.-,_(( ( ( (  (           _         )      ). . \  ',  '!    \^       |                           |   | |
' |  '", ( (  (  (          ', . ..:/     .::     ;   '",'!, __)^      |   .                       |   | |
| |     ',_ (  (    .         :.:.:/    .::  :"    \     ','',,,'      |    .                      |   | |
| |        ""-,_ .(.;         |:::/  _.:   +` '':   ;       '"""       |     .             5σ      |   | |
| |             ",_.:\        |::(  (:  :``      '-._'=-,              |      .                    |   | |
| |                ',_)       ::::',_', ",           "\)\)             |       '.   ..             |   | |
| |                  /        ;___.`)\\)  )                            |         '.'  '            |   | |
| |                 `         `    /     /                             |               '.          |   | |
| |                :        .`   .`    .`                              |                 '......   |   | |
| |                : . . . :    /. . ./                                |___________________________|   | |
| |                :. . ..`    :. . ./                                                                 | |
| |                /...:/     |:.:..`   _________________              ______  _  _    _               | |
| |               ;:::::       ::::(   |_____   ___   __ \            |  ____|(_)| |  | |              | |
| |                :::::        )::::,__     | |   | |__) | ___ __  __| |__    _ | |_ | |_  ___  _ __  | |
| |                ::::;        :::::::::;,  | |   |  _  / / _ \\ \/ /|  __|  | || __|| __|/ _ \| '__| | |
| |               /:::::'--=,  /`.-""--,;    | |   | | \ \|  __/ >  < | |     | || |_ | |_|  __/| |    | |
| |              /::::::::."-\ |"       "\   |_|   |_|  \_\\___|/_/\_\|_|     |_| \__| \__|\___||_|    | |
| |_____________(::(`'"-,;_____________________________________________________________________________| |
|________________|"______"\______________________________________________________________________________|
EOF

echo -e "\n"${name}"\nversion: "${version}"\n\nrun name: "${run_name}"\n"

timestamp_run="$(date -u "+%Y-%m-%dT%H%M%S")Z"
filename_configuration="config/"${run_name}".config"
filename_histograms_log=""${timestamp_run}"_histograms_log.txt"
filename_prefit_plots_log=""${timestamp_run}"_prefit_plots_log.txt"
filename_workspace_log=""${timestamp_run}"_workspace_log.txt"
filename_fit_log=""${timestamp_run}"_fit_log.txt"
filename_limit_log=""${timestamp_run}"_limit_log.txt"
filename_significance_log=""${timestamp_run}"_significance_log.txt"

if [ ! -f "${filename_configuration}" ]; then
    echo -e "error: configuration "${filename_configuration}" not found\n"
    exit 1
fi

echo "run start $(date -u "+%Y-%m-%dT%H%M%S")Z"

rm -rf "${run_name}"

echo "access input histograms $(date -u "+%Y-%m-%dT%H%M%S")Z"
time ./myFit.exe h "${filename_configuration}" > >(tee "${filename_histograms_log}")

echo "draw pre-fit plots $(date -u "+%Y-%m-%dT%H%M%S")Z"
time ./myFit.exe d "${filename_configuration}" > >(tee "${filename_prefit_plots_log}")

echo "create the RooStats XMLs and workspace $(date -u "+%Y-%m-%dT%H%M%S")Z"
time ./myFit.exe w "${filename_configuration}" > >(tee "${filename_workspace_log}")

echo "fit the workspace $(date -u "+%Y-%m-%dT%H%M%S")Z"
time ./myFit.exe f "${filename_configuration}" > >(tee "${filename_fit_log}")

echo "exclusion limit $(date -u "+%Y-%m-%dT%H%M%S")Z"
time ./myFit.exe l "${filename_configuration}" > >(tee "${filename_limit_log}")

echo "significance $(date -u "+%Y-%m-%dT%H%M%S")Z"
time ./myFit.exe s "${filename_configuration}" > >(tee "${filename_significance_log}")

echo "run stop $(date -u "+%Y-%m-%dT%H%M%S")Z"
