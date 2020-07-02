#!/bin/bash
cd $1
sh Make.sh style
mv style_*.h $2
