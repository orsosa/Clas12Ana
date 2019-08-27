#!/bin/bash
if [ $# -ne 1 ]; then
    echo "one parameter name must be given."
    exit 1
fi
echo "looking "$1" on varList"
grep "varList =" examples/get_simple_tuple_th_mf.cxx | sed 's/:/\n/g' | grep -n $1

echo "looking "$1" on varListElec"
grep "varListElec =" examples/get_simple_tuple_th_mf.cxx | sed 's/:/\n/g' | grep -n $1
