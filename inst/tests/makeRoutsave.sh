#!/bin/bash

rm -f *.Rout.save
for i in *.R
do
  R --no-save < $i > $i'out.save' 2>&1
done

