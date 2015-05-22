# this script is an example for turning the tables from the fitter into latex format

while read line
do
  line=${line/'|'/' $'}
  line=${line//'|'/'$ & $'}
  line=${line//'+/-'/'\pm'}
  line=${line//'#'/'\'}
  line="${line%?}"
  line="${line%?}"
  line="${line%?}"
  line+=' \\'
  echo $line
done < $1
