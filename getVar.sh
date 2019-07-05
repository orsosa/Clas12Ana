#!/bin/bash
if [ $# -ne 1 ]; then
    echo "one parameter name must be given."
    exit 1
fi
echo "looking "$1
grep "VarList =" examples/get_simple_tuple.cxx | sed 's/:/\n/g' | grep -n $1
