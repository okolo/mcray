import numpy as np
import sys

if len(sys.argv)!=3:
    print('usage: python', sys.argv[0],'tr.log.sorted recordId', file=sys.stderr)
    print('tr.log.sorted is produced by command sort -g -r -k 4 tr.log > tr.log.sorted', file=sys.stderr)
    sys.exit(1)

if len(sys.argv)>2:
    id = int(sys.argv[2])

data = np.loadtxt(sys.argv[1])


indexes = data[:, 0] == id
subrecords = data[indexes]
if len(subrecords) == 0:
    print("id",id,'not found',file=sys.stderr)
    sys.exit(1)
par_id = subrecords[0,1]
assert np.sum(subrecords[:,1]!=par_id) == 0 # make sure all records have the same parent

records = [subrecords]

while par_id>0:
    t_min = np.min(subrecords[:, 3])
    indexes = np.logical_and(data[:, 0] == par_id, data[:,3] <= t_min)
    subrecords = data[indexes]
    if len(subrecords)==0:
        break
    par_id = subrecords[0, 1]
    assert np.sum(subrecords[:,1]!=par_id) == 0 # make sure all records have the same parent

    records += [subrecords]
    #data = data[np.logical_not(indexes)]

r = np.vstack(records)[::-1,:]
with open(str(int(id)),'wt') as out:
    for l in r:
        print(int(l[0]),int(l[1]),'\t'.join(["{0:g}".format(n) for n in l[2:]]), file=out)

