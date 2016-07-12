#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
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

usage:
    run.py [options]

options:
    -h, --help           display help message
    --version            display version and exit
    --runname=NAME       run name                 [default: ttHbb]
    --histograms=BOOL    run histograms           [default: true]
    --prefitplots=BOOL   run prefit plots         [default: true]
    --workspace=BOOL     run workspace            [default: true]
    --fit=BOOL           run fit                  [default: true]
    --limit=BOOL         run limit                [default: true]
    --significance=BOOL  run significance         [default: false]
"""

name    = "TRexFitter run"
version = "2016-07-12T1217Z"
logo =\
'''
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
| |                  /        ;___.`)\\\)  )                            |         '.'  '            |   | |
| |                 `         `    /     /                             |               '.          |   | |
| |                :        .`   .`    .`                              |                 '......   |   | |
| |                : . . . :    /. . ./                                |___________________________|   | |
| |                :. . ..`    :. . ./                                                                 | |
| |                /...:/     |:.:..`   _________________              ______  _  _    _               | |
| |               ;:::::       ::::(   |_____   ___   __ \            |  ____|(_)| |  | |              | |
| |                :::::        )::::,__     | |   | |__) | ___ __  __| |__    _ | |_ | |_  ___  _ __  | |
| |                ::::;        :::::::::;,  | |   |  _  / / _ \\\ \/ /|  __|  | || __|| __|/ _ \| '__| | |
| |               /:::::'--=,  /`.-""--,;    | |   | | \ \|  __/ >  < | |     | || |_ | |_|  __/| |    | |
| |              /::::::::."-\ |"       "\   |_|   |_|  \_\\\___|/_/\_\|_|     |_| \__| \__|\___||_|    | |
| |_____________(::(`'"-,;_____________________________________________________________________________| |
|________________|"______"\______________________________________________________________________________|
'''

import datetime
import docopt
import os
import subprocess
import sys
import time

def main(options):

    run_name         = options["--runname"]
    run_histograms   = string_to_bool(options["--histograms"])
    run_prefit_plots = string_to_bool(options["--prefitplots"])
    run_workspace    = string_to_bool(options["--workspace"])
    run_fit          = string_to_bool(options["--fit"])
    run_limit        = string_to_bool(options["--limit"])
    run_significance = string_to_bool(options["--significance"])

    print(logo)

    print("\n{name}\nversion: {version}\n\nrun name: {run_name}\n".format(
        name     = name,
        version  = version,
        run_name = run_name
    ))

    print("run options:")
    print("""
- run histograms:   {run_histograms}
- run prefit plots: {run_prefit_plots}
- run workspace:    {run_workspace}
- run fit:          {run_fit}
- run limit:        {run_limit}
- run significance: {run_significance}""".format(
        run_histograms   = run_histograms,
        run_prefit_plots = run_prefit_plots,
        run_workspace    = run_workspace,
        run_fit          = run_fit,
        run_limit        = run_limit,
        run_significance = run_significance
    ))

    datetime_object_current_time_UTC = datetime.datetime.utcnow()

    timestamp_run             = datetime_object_current_time_UTC.strftime("%Y-%m-%dT%H%M%SZ")
    filename_configuration    = "config/{run_name}.config".format(run_name = run_name)
    filename_histograms_log   = "{timestamp_run}_histograms_log.txt".format(timestamp_run = timestamp_run)
    filename_prefit_plots_log = "{timestamp_run}_prefit_plots_log.txt".format(timestamp_run = timestamp_run)
    filename_workspace_log    = "{timestamp_run}_workspace_log.txt".format(timestamp_run = timestamp_run)
    filename_fit_log          = "{timestamp_run}_fit_log.txt".format(timestamp_run = timestamp_run)
    filename_limit_log        = "{timestamp_run}_limit_log.txt".format(timestamp_run = timestamp_run)
    filename_significance_log = "{timestamp_run}_significance_log.txt".format(timestamp_run = timestamp_run)

    executable                = "./myFit.exe"

    if not os.path.isfile(filename_configuration):
        print("error: configuration {filename} not found\n".format(
            filename = filename_configuration
        ))
        sys.exit()

    if not os.path.isfile(executable):
        print("error: executable {filename} not found\n".format(
            filename = executable
        ))
        sys.exit()

    print("run start {timestamp}".format(
        timestamp = timestamp_run
    ))

    print(
        """
log files:

- {filename_histograms_log}
- {filename_prefit_plots_log}
- {filename_workspace_log}
- {filename_fit_log}
- {filename_limit_log}
- {filename_significance_log}\n""".format(
        filename_histograms_log   = filename_histograms_log,
        filename_prefit_plots_log = filename_prefit_plots_log,
        filename_workspace_log    = filename_workspace_log,
        filename_fit_log          = filename_fit_log,
        filename_limit_log        = filename_limit_log,
        filename_significance_log = filename_significance_log
    ))

    if os.path.isdir(run_name):
        directory_name_tmp = "{timestamp}_backup_{run_name}".format(
            timestamp = timestamp_run,
            run_name  = run_name
        )
        print("existing results found at directory {directory_name} -- move to directory {directory_name_tmp}\n".format(
            directory_name     = run_name,
            directory_name_tmp = directory_name_tmp
        ))
        command = "mv {directory_name} {directory_name_tmp}".format(
            directory_name     = run_name,
            directory_name_tmp = directory_name_tmp
        )
        engage_command(command)

    command = ""
    if run_histograms == True:
        command = command + """
echo "access input histograms $(date -u "+%Y-%m-%dT%H%M%S")Z"
time {executable} h {filename_configuration} > >(tee {filename_histograms_log})
"""
    if run_prefit_plots == True:
        command = command + """
echo "draw pre-fit plots $(date -u "+%Y-%m-%dT%H%M%S")Z"
time {executable} d {filename_configuration} > >(tee {filename_prefit_plots_log})
"""
    if run_workspace == True:
        command = command + """
echo "create the RooStats XMLs and workspace $(date -u "+%Y-%m-%dT%H%M%S")Z"
time {executable} w {filename_configuration} > >(tee {filename_workspace_log})
"""
    if run_fit == True:
        command = command + """
echo "fit the workspace $(date -u "+%Y-%m-%dT%H%M%S")Z"
time {executable} f {filename_configuration} > >(tee {filename_fit_log})
"""
    if run_limit == True:
        command = command + """
echo "exclusion limit $(date -u "+%Y-%m-%dT%H%M%S")Z"
time {executable} l {filename_configuration} > >(tee {filename_limit_log})
"""
    if run_significance == True:
        command = command + """
echo "significance $(date -u "+%Y-%m-%dT%H%M%S")Z"
time {executable} s {filename_configuration} > >(tee {filename_significance_log})
"""
    command = command + """
echo "run stop $(date -u "+%Y-%m-%dT%H%M%S")Z"
echo "run time: {timestamp_run}--$(date -u "+%Y-%m-%dT%H%M%S")Z"
"""
    command = command.format(
        timestamp_run             = timestamp_run,
        executable                = executable,
        filename_configuration    = filename_configuration,
        filename_histograms_log   = filename_histograms_log,
        filename_prefit_plots_log = filename_prefit_plots_log,
        filename_workspace_log    = filename_workspace_log,
        filename_fit_log          = filename_fit_log,
        filename_limit_log        = filename_limit_log,
        filename_significance_log = filename_significance_log
    )

    print_line(character = ".")
    print("execute command:\n{command}".format(
        command = command
    ))
    print_line(character = ".")

    print("\nDo you want to continue with the run? [y/n]")
    response = get_y_or_n()
    if response == "n":
        print("\nexit")
        print_line()
        sys.exit()
    else:
        print("\nstart run\n")

    engage_command(command)

    print_line()

def get_keystroke():
    import sys
    import tty
    import termios
    fd = sys.stdin.fileno()
    old_settings = termios.tcgetattr(fd)
    try:
        tty.setraw(sys.stdin.fileno())
        character = sys.stdin.read(1)
    finally:
        termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
    return character

def get_y_or_n():
    character = None
    while character != "y" and character != "n":
        character = get_keystroke().lower()
    return character

def interrogate(
    prompt   = None,
    default  = None,
    template = "{prompt} [default: {default}]: "
    ):
    message = template.format(
        prompt  = prompt,
        default = default,
    )
    response = get_input(message)
    if response:
        return response
    else:
        return default

def pause(
    prompt = "\nPress Enter to continue."
    ):
    get_input(prompt)

def get_input(
    prompt = None
    ):
    if sys.version_info >= (3, 0):
        return input(prompt)
    else:
        return raw_input(prompt)

def engage_command(
    command = None
    ):
    process = subprocess.Popen(
        [command],
        shell      = True,
        executable = "/bin/bash")
    process.wait()
    output, errors = process.communicate()
    return output

def terminal_width():
    return(
        int(
            subprocess.Popen(
                ["tput", "cols"],
                stdout = subprocess.PIPE
            ).communicate()[0].decode("utf-8").strip("\n")
        )
    )

def line_string(
    character = "-"
    ):
    _terminal_width = terminal_width()
    line = ""
    for column in range(0, _terminal_width):
        line += character
    return(line)

def print_line(
    character = "-"
    ):
    print(line_string(
        character = character
    ))

def string_to_bool(x):
    return x.lower() in ("yes", "true", "t", "1")

if __name__ == "__main__":
    options = docopt.docopt(__doc__)
    if options["--version"]:
        print(version)
        exit()
    main(options)
