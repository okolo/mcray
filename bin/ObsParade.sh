N=300000
THREADS=62

echo "Test l EGMF 1"
./CRbeam -p 10 -z 0.00245 -E 1e19 -e 1e15 -mft -mf 1e-9 -nth $THREADS -oo -o output/articleVar_at_z_0.0024_3 -jet 0.03 -N $N -rseed 42 --lEGMF 1 -obs 1 -obsr 100 -os -jett 0.1343084095931479 -jetp 6.216272108729293 -obsz 10000
echo "Test different positions"
./CRbeam -p 10 -z 0.00234 -E 1e19 -e 1e15 -mft -mf 1e-9 -nth $THREADS -oo -o output/articleVar_at_D9000 -jet 0.2 -N $N -rseed 42 --lEGMF 1 -obs 1 -obsr 100 -os -jett 0.1 -jetp 6.2 -obsz 9000
./CRbeam -p 10 -z 0.00234 -E 1e19 -e 1e15 -mft -mf 1e-9 -nth $THREADS -oo -o output/articleVar_at_D8000 -jet 0.2 -N $N -rseed 42 --lEGMF 1 -obs 1 -obsr 100 -os -jett 0.1 -jetp 6.2 -obsz 8000
./CRbeam -p 10 -z 0.00234 -E 1e19 -e 1e15 -mft -mf 1e-9 -nth $THREADS -oo -o output/articleVar_at_D7000 -jet 0.2 -N $N -rseed 42 --lEGMF 1 -obs 1 -obsr 100 -os -jett 0.1 -jetp 6.2 -obsz 7000
./CRbeam -p 10 -z 0.00234 -E 1e19 -e 1e15 -mft -mf 1e-9 -nth $THREADS -oo -o output/articleVar_at_D5000 -jet 0.3 -N $N -rseed 42 --lEGMF 1 -obs 1 -obsr 100 -os -jett 0.1 -jetp 6.2 -obsz 5000
./CRbeam -p 10 -z 0.00234 -E 1e19 -e 1e15 -mft -mf 1e-9 -nth $THREADS -oo -o output/articleVar_at_D1000 -jet 0.4 -N $N -rseed 42 --lEGMF 1 -obs 1 -obsr 100 -os -jett 0.1 -jetp 6.2 -obsz 1000
./CRbeam -p 10 -z 0.00234 -E 1e19 -e 1e15 -mft -mf 1e-9 -nth $THREADS -oo -o output/articleVar_at_D800 -jet 0.8 -N $N -rseed 42 --lEGMF 1 -obs 1 -obsr 100 -os -jett 0.1 -jetp 6.2 -obsz 800
#./CRbeam -p 10 -z 0.00234 -E 1e19 -e 1e15 -mft -mf 1e-9 -nth $THREADS -o output/articleVar_at_D500 -jet 0.9 -N 90000 -rseed 42 --lEGMF 1 -obs 1 -obsr 100 -os -jett 0.1 -jetp 6.2 -obsz 500
#./CRbeam -p 10 -z 0.00234 -E 1e19 -e 1e15 -mft -mf 1e-9 -nth 48 -o output/articleVar_at_D200 -jet 1.1 -N 90000 -rseed 42 --lEGMF 1 -obs 1 -obsr 100 -os -jett 0.1 -jetp 6.2 -obsz 200
./CRbeam -p 10 -z 0.00234 -E 1e19 -e 1e15 -mft -mf 1e-9 -nth $THREADS -oo -o output/articleVar_at_D100 -jet 2.2 -N 90000 -rseed 42 --lEGMF 1 -obs 1 -obsr 100 -os -jett 0.1 -jetp 6.2 -obsz 100
