#!/bin/bash

while true;
do
	d=$(date)
	#r=$(ping -c 60 ya.ru | grep "packet loss")
	r=$(ping -c 60 172.20.91.18 | grep "packet loss")
	echo $d";"$r
done
