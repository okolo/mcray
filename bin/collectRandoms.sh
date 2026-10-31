#scp -r "cluster83:output/srcRand*" /tmp/dolgikh/
scp -r "cluster84:output/srcRand37_*" /tmp/dolgikh37/
#scp -r "cluster82:output/srcRand*" /tmp/dolgikh/
scp -r "cluster81:output/srcRand37_*" /tmp/dolgikh37/
#scp -r "cluster54:output/srcRand*" /tmp/dolgikh/
#scp -r "cluster53:output/srcRand*" /tmp/dolgikh/
cp -r output/srcRand37_* /tmp/dolgikh37/
for var in $(ls /tmp/dolgikh37/)
do
	selectInTarget /tmp/dolgikh37/$var/z0/proton > /tmp/dolgikh37/$var/z0/proton_target
#	generateHPmap.py /tmp/dolgikh/$var/z0/proton_target
done
