#!/bin/bash

echo "########################################"
echo "pydeface"
echo "$1"
echo "$2"
echo "########################################"

pydeface "$1" --outfile "$2"
