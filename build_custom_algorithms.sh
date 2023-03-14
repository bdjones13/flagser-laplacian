#!/bin/sh
cd `git rev-parse --show-cdup`
if [ 0 -lt `git status -s algorithms.math | wc -l` ]
then
    echo "You modified algorithms.math so I install the required python packages to compile them. If you want to prevent this, execute git checkout algorithms.math."
    python3 -m pip install --user -r requirements.txt
    python3 ./parser/math_parser.py
fi
# commented by Ben Jones 2023-03-01
