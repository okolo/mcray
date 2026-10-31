#!/bin/bash
echo "13">13
./SOPHIA_CS <13 | awk 'NR>3' >p
echo "14">14
./SOPHIA_CS <14 | awk 'NR>3' >n
rm 13 14
