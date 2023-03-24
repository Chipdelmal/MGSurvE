#!/bin/bash

for m in "sum" "man" "max"; do
    echo "* Optimizing $m"
    for i in {1..30}; do
        echo "* Running iteration $i"
        python YKN-Discrete.py "YKN" "$m" "$i"
    done
done
