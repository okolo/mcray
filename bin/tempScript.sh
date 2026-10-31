#!/usr/bin/bash
echo "Sleeping..."
sleep 43200
echo "Killing CRbeam"
kill 127936
sleep 10
echo "Kill main script"
kill 29175
sleep 10
echo "Call Flux analyzer"
cd /net/minus/home/dolgikh/SemikozPO/CRbeam_v0/bin
python FluxBatchGetter.py
