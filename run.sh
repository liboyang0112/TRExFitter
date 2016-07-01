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
version="2016-06-30T1444Z"

if [ -z "${1}" ]; then
    run_name="ttHbb"
else
    run_name="${1}"
fi

main(){
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

executable="./myFit.exe"

if [ ! -f "${filename_configuration}" ]; then
    echo -e "error: configuration "${filename_configuration}" not found\n"
    exit 1
fi

if [ ! -f "${executable}" ]; then
    echo -e "error: executable "${executable}" not found\n"
    exit 1
fi

print_line

echo "run start $(date -u "+%Y-%m-%dT%H%M%S")Z"

echo -e "\nlog files:\n"
echo -e "- "${filename_histograms_log}""
echo -e "- "${filename_prefit_plots_log}""
echo -e "- "${filename_workspace_log}""
echo -e "- "${filename_fit_log}""
echo -e "- "${filename_limit_log}""
echo -e "- "${filename_significance_log}"\n"

if [ -d "${run_name}" ]; then
    tmp_directory=""$(date -u "+%Y-%m-%dT%H%M%S")Z"_backup_"${run_name}""
    echo -e "existing results found at directory "${run_name}" -- move to directory "${tmp_directory}"\n"
    mv "${run_name}" "${tmp_directory}"
fi

echo "access input histograms $(date -u "+%Y-%m-%dT%H%M%S")Z"
time "${executable}" h "${filename_configuration}" > >(tee "${filename_histograms_log}")

echo "draw pre-fit plots $(date -u "+%Y-%m-%dT%H%M%S")Z"
time "${executable}" d "${filename_configuration}" > >(tee "${filename_prefit_plots_log}")

echo "create the RooStats XMLs and workspace $(date -u "+%Y-%m-%dT%H%M%S")Z"
time "${executable}" w "${filename_configuration}" > >(tee "${filename_workspace_log}")

echo "fit the workspace $(date -u "+%Y-%m-%dT%H%M%S")Z"
time "${executable}" f "${filename_configuration}" > >(tee "${filename_fit_log}")

echo "exclusion limit $(date -u "+%Y-%m-%dT%H%M%S")Z"
time "${executable}" l "${filename_configuration}" > >(tee "${filename_limit_log}")

echo "significance $(date -u "+%Y-%m-%dT%H%M%S")Z"
time "${executable}" s "${filename_configuration}" > >(tee "${filename_significance_log}")

echo "run stop $(date -u "+%Y-%m-%dT%H%M%S")Z"

echo "run time: "${timestamp_run}"--$(date -u "+%Y-%m-%dT%H%M%S")Z"

print_line
}

#¯`·.¸¸.·´¯`·.¸¸.·´¯`·.¸¸.·´¯`·.¸¸.·´¯`·.¸¸.·´¯`·.¸¸.·´¯`·.¸¸.·´¯`·.¸¸.·´><(((º>
print_line(){
################################################################################
if [ "${1}" = "--interrogate" ]; then
IFS= read -d '' function_information << "EOF"
This function prints one line on the terminal. It may be used for the purpose of
terminal output legibility. This function initially determines the terminal
width and then prints one line of _ characters.
EOF
return
fi
################################################################################
    number_of_lines=1
    terminal_width="$(return_terminal_dimension "width")"
    number_of_characters_to_print=$(\
        echo "${number_of_lines}*(${terminal_width})-1" | bc\
    )
    print_character "_" "${number_of_characters_to_print}"
}
 
#¯`·.¸¸.·´¯`·.¸¸.·´¯`·.¸¸.·´¯`·.¸¸.·´¯`·.¸¸.·´¯`·.¸¸.·´¯`·.¸¸.·´¯`·.¸¸.·´><(((º>
print_character(){
################################################################################
if [ "${1}" = "--interrogate" ]; then
IFS= read -d '' function_information << "EOF"
This function takes two arguments, the first being the character to print and
the second being the number of times to print the character. This function
prints the specified character a specified number of times without carriage
returns.
EOF
return
fi
################################################################################
    character="${1}"
    number_of_times_to_print_character="${2}"
    for (( \
        current_print_number = 0; \
        current_print_number<=${number_of_times_to_print_character}; \
        current_print_number++ \
        )); do
        echo -n "${character}"
    done
}
 
#¯`·.¸¸.·´¯`·.¸¸.·´¯`·.¸¸.·´¯`·.¸¸.·´¯`·.¸¸.·´¯`·.¸¸.·´¯`·.¸¸.·´¯`·.¸¸.·´><(((º>
return_terminal_dimension(){
################################################################################
if [ "${1}" = "--interrogate" ]; then
IFS= read -d '' function_information << "EOF"
This function takes one argument, the required dimension of the
terminal. The required dimension then is returned. If no argument is
specified, nothing is returned. The possible dimensions are as follows:
- size:   the width and height of the terminal separated by a space (e.g. 24 80)
- height: the height of the terminal (e.g. 24)
- width:  the width of the terminal (e.g. 80)
EOF
return
fi
################################################################################
    dimension="${1}"
    if [ "${dimension}" = "size" ]; then
        stty size
    else
        if [ "${dimension}" = "height" ]; then
            stty size | cut -d" " -f1
        else
            if [ "${dimension}" = "width" ]; then
                stty size | cut -d" " -f2
            fi
        fi
    fi
}

main