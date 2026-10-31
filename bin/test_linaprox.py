import numpy as np
import matplotlib.pyplot as plt

print("Load MFs")
interpolated = np.loadtxt("output/MFsaveTest_TestInterpolation2/magneticField.out")
lowquality = np.loadtxt("output/MFsaveTest_orig_seed42/magneticField.out")

print(interpolated.shape[0]/lowquality.shape[0])

def plot_slice(arr,x=None,y=None,z=None,show=True,title=""):
	sel = np.ones_like(arr[:,0]).astype(bool)
	x_axis_ind = 0
	if x is not None:
		sel = sel & (arr[:,0]==x)
	else:
		x_axis_ind = 0
	if y is not None:
		sel = sel & (arr[:,1]==y)
	else:
		x_axis_ind = 1
	if z is not None:
		sel = sel & (arr[:,2]==z)
	else:
		x_axis_ind = 2
	x_axis_labels = ["x, Mpc","y, Mpc", "z, Mpc"]
	mslice = arr[ sel ]
	plt.plot(mslice[:,x_axis_ind],mslice[:,3]*1e9,'.-',label="Bx, nG")
	plt.plot(mslice[:,x_axis_ind],mslice[:,4]*1e9,'.-',label="By, nG")
	plt.plot(mslice[:,x_axis_ind],mslice[:,5]*1e9,'.-',label="Bz, nG")
	plt.title(f"x={x} Mpc, y={y} Mpc, z={z} Mpc, {title}")
	plt.legend()
	if(show):plt.show()

interpolated_coords = list(set(interpolated[:,0]))
lowquality_coords = list(set(lowquality[:,0]))
interpolated_coords.sort()
lowquality_coords.sort()

# BORDER, Nearest approx

plot_slice(interpolated,x=interpolated_coords[0],y=interpolated_coords[0],show=False)
plot_slice(lowquality,x=lowquality_coords[0],y=lowquality_coords[0],title="border")

plot_slice(interpolated,y=interpolated_coords[0],z=interpolated_coords[0],show=False)
plot_slice(lowquality,y=lowquality_coords[0],z=lowquality_coords[0],title="border")

plot_slice(interpolated,x=interpolated_coords[0],z=interpolated_coords[0],show=False)
plot_slice(lowquality,x=lowquality_coords[0],z=lowquality_coords[0],title="border")

# inner, linear approx
interpolated_linear_test = interpolated_coords[5]
ind = abs(np.array(lowquality_coords)-interpolated_linear_test).argmin()

print("interpolated_linear_test =",interpolated_linear_test)
print("lowquality_coords[",ind,"] =",lowquality_coords[ind])

plot_slice(interpolated,x=interpolated_coords[5],y=interpolated_coords[5],show=False)
plot_slice(lowquality,x=lowquality_coords[ind],y=lowquality_coords[ind])

plot_slice(interpolated,y=interpolated_coords[5],z=interpolated_coords[5],show=False)
plot_slice(lowquality,y=lowquality_coords[ind],z=lowquality_coords[ind])

plot_slice(interpolated,x=interpolated_coords[5],z=interpolated_coords[5],show=False)
plot_slice(lowquality,x=lowquality_coords[ind],z=lowquality_coords[ind])
