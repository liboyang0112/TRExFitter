#!/bin/bash
# this script is an example for turning the tables from the fitter into latex format

##____________________________________________
## Converts all lines into latex format
##    argument 1 is the input table
##    argument 2 is the output table
function convert_all_lines {
    nline=0
    total_lines=`wc $1 | awk '{print $1}'`
    isSyst=`echo $1 | grep -b -o syst`
    if [[ "$isSyst" = "" ]]; then
        isSyst=0
    else
        isSyst=1
    fi

    while read line
        do
            if [[ $isSyst -eq 0 ]]; then
                if [ $nline -eq 0 ] || [ $nline -eq $((total_lines-1)) ] || [ $nline -eq $((total_lines-2)) ] ; then
                    echo "\t \t \\hline" >> $2
                fi
            else
                if [ $nline -eq 0 ] ; then
                    echo "\t \t \\hline" >> $2
                fi
            fi
            line=${line/'|'/' $'}
            line=${line//'|'/'$ & $'}
            line=${line//'pm'/'\\pm'}
            line=${line//'#'/'\'}
            line="\t \t ${line%?}"
            line="${line%?}"
            line="${line%?}"
            line+=' \\\\'
            echo $line >> $2
            nline=$((nline+1))
            if [ $nline -eq 1 ] || [ $nline -eq $total_lines ] ; then
                echo "\t \t \\hline" >> $2
            fi
        done < $1
}

##____________________________________________
## Write the preample of the latex tables
function write_preamble {
    column=0
    while read line
        do
        column=`grep -o "[|]" <<<"$line" | wc -l`
        column=$((column-1))
        break
        done < $1
    echo "\\\begin{table}" > $2
    line="\t \\\begin{tabular}{| l |"
    for i in $(seq 1 $column)
        do
            line+=" c |"
        done
    line+="}"
    echo $line >> $2
}

##____________________________________________
## Write the end of the latex tables
function write_end {
    echo "\t \\\end{tabular}" >> $2
    echo "\\\end{table}" >> $2
}

##____________________________________________
## Build the latex tables (assuming one provided a config files)
function tableToLatex {
    if [[ $# -eq 0 ]]; then
            echo "=> No config file provided. Plase provide one."
    else
            fitName=`less $1 | grep "Job:" | awk '{print $2}'`
            fitName=${fitName//'"'/''}
            for inFile in `find $fitName/Tables -name '*.txt'`
                do
                    echo "-> Dumping info from $inFile"
                    outFile=${inFile/".txt"/".tex"}
                    write_preamble $inFile $outFile
                    convert_all_lines $inFile $outFile
                    write_end $inFile $outFile
                done
    fi
}

tableToLatex $1