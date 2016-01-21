#!/bin/bash
COLUMNS=$(tput cols)
export COLUMNS

java -jar $(ls dist/*.jar) ${@} 

