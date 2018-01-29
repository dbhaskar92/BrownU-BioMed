#!/bin/bash
# Author: Dhananjay Bhaskar

# Delete CSV file
rm -rf *.csv

# Delete all subdirectories
find . -maxdepth 1 -mindepth 1 -type d -exec rm -rf {} \;
