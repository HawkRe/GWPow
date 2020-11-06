import h5py
import numpy as np
from scipy import fftpack
import sys
from fractions import Fraction

lattice_dim = 192
boxsize = 1600000.0
delta_interval =  boxsize/lattice_dim
K_IR = 2*np.pi/boxsize 
coefficient = np.power(K_IR,3)*np.power(delta_interval,6)/(np.power(4*np.pi,4)*np.power(boxsize,3))

def XYZ():
    '''
    real x,y,z pattern
    '''
    points = np.arange(0,192*192*192)
    #quarter = 48
    index = 3520415 # center + center*N + center*N**2
    points[index],points[-1] = 1778736,1778736 # to solve the problem of dividing by zero
    return points.reshape(192,192,192)

# faster effective momentum pattern
def Lambda(value,N,i,j,l,m):
    def projector(value,N,i,j):
        def delta_function(i,j):
            '''
            \delta_{i,j} function  
            '''
            if (i == j):
                return 1
            else: return 0
        def eff_momentum(x):
            eff = np.sin(2*np.pi*x/N)
            return eff
        def coordinate(value,N):
            #! Seems Problematic, but compatible with row-major fftw-3 and numpy conventions!
            center = int((N-1)/2)
            x = value//(N**2) 
            y = (value - x*N**2)//N 
            z = value - y*N - x*N**2 # there is a issue: the divider might be zero! so we must deal with 3-d np-array ahead-of-time
            return eff_momentum(x-center),eff_momentum(y-center),eff_momentum(z-center)
        def norm(vector):
            return vector[0]**2 + vector[1]**2 + vector[2]**2
        projector = delta_function(i-1,j-1) - coordinate(value,N)[i-1]*coordinate(value,N)[j-1]/norm(coordinate(value,N))
        return projector

    eff_lambda = projector(value,N,i,l)*projector(value,N,j,m) - 0.5*projector(value,N,i,j)*projector(value,N,l,m) 
    # we need to add extra: -0.5(\delta_{il}+\delta_{jm}) + 0.25(\delta_{ij}+\delta_{lm}) + 0.125 with all other elements being zero except [95][95][95], [191][191][191]
    return eff_lambda

# We assume that the data of perturbation has been distracted from the hdf5 file.
#! Recommend to use numpy-array intrinsic multiplication
def sum_modes(value,N,dataset):
    def Extra(N,i,j,l,m):
        def delta_function(i,j):
            '''
            \delta_{i,j} function  
            '''
            if (i == j):
                return 1
            else: return 0
            
        eff_value = -(delta_function(i,l)+delta_function(j,m))/6 + (delta_function(i,j)+delta_function(l,m))/12 + 5/72
        center = int((N-1)/2)
        index = center + center*N + center*N**2
        EXTRA = np.zeros(N*N*N) 
        EXTRA[index],EXTRA[-1] = eff_value, eff_value 
        return EXTRA.reshape(N,N,N)

    def gauge_component(value,i,j,N,dataset):
        eff = (Lambda(value,N,i,j,1,1)+Extra(N,i,j,1,1))*dataset['w11'] + (Lambda(value,N,i,j,1,2)+Extra(N,i,j,1,2))*dataset['w12'] \
            + (Lambda(value,N,i,j,1,3)+Extra(N,i,j,1,3))*dataset['w13'] + (Lambda(value,N,i,j,2,2)+Extra(N,i,j,2,2))*dataset['w22'] \
            + (Lambda(value,N,i,j,2,3)+Extra(N,i,j,2,3))*dataset['w23'] + (Lambda(value,N,i,j,3,3)+Extra(N,i,j,3,3))*dataset['w33']
        return eff
    sum = np.abs(np.power(gauge_component(value,1,1,N,dataset),2)) + np.abs(np.power(gauge_component(value,1,2,N,dataset),2)) + np.abs(np.power(gauge_component(value,1,3,N,dataset),2)) \
        + np.abs(np.power(gauge_component(value,2,2,N,dataset),2)) + np.abs(np.power(gauge_component(value,2,3,N,dataset),2)) + np.abs(np.power(gauge_component(value,3,3,N,dataset),2))
    return sum

def gw_fft(dataset):
    '''
    We require the dataset to be dictionary object: the key is the name of variables (w11 ... w33), 
    and the value is the data 
    we have shift the zero frequency to the "center" position
    '''
    dataset['w11'] = fftpack.fftshift(fftpack.fftn(dataset['w11'].astype('float64')))
    dataset['w12'] = fftpack.fftshift(fftpack.fftn(dataset['w12'].astype('float64')))
    dataset['w13'] = fftpack.fftshift(fftpack.fftn(dataset['w13'].astype('float64')))
    dataset['w22'] = fftpack.fftshift(fftpack.fftn(dataset['w22'].astype('float64')))
    dataset['w23'] = fftpack.fftshift(fftpack.fftn(dataset['w23'].astype('float64')))
    dataset['w33'] = fftpack.fftshift(fftpack.fftn(dataset['w33'].astype('float64')))
    return dataset


def gw_spectrum(sum_power,N):
    '''
    calculate the power-spectrum of gravitational waves, input dataset to be dictionary-like object.
    '''
    grid_h = int((N-1)/2)
    distance = np.sqrt(3) # describing the searching meta-length
    max_k = int(N/(distance*2)) # the logic is that: we'll later shift the frequency distribution with the zero frequency to the center
    power_spectrum = []
    eff_spectrum = 0.0
    for i in range(1,max_k):
        norm = 0 
        for x_1 in range(-int(i*distance),int(i*distance)+1):
            for y_1 in range(-int(i*distance),int(i*distance)+1):
                for z_1 in range(-int(i*distance),int(i*distance)+1):
                    if (3*(i-1)*(i-1) <= (x_1**2 + y_1**2 + z_1**2) < 3*i*i):
                        norm = norm + sum_power[grid_h + x_1][grid_h + y_1][grid_h + z_1]
                    else:
                        pass
        eff_spectrum = norm * np.power((i-0.5)*3,0.5)
        power_spectrum.append(eff_spectrum)
    return power_spectrum

def manipulate_data(size,time,step_interval,main):
    def treat_string(string_step):
        '''
        be compatible with SAMRAI output format of folder and files, string_step is kind of string type
        '''
        # more faster reverse string = string[::-1]
        real_string = '0'*(5-len(string_step)) + string_step 
        return real_string 
    def extract_data(directory,grid_size):
        summary_file = directory + '/' + "summary.samrai"
        basic_info = h5py.File(summary_file,'r')
        num_patches = int(basic_info["/BASIC_INFO/number_patches_at_level"][...])
        patch_extents = basic_info["/extents/patch_extents"][...]
        vars = list(np.char.decode(basic_info["/BASIC_INFO/var_names"][...]))
        data = {}
        for j in range(0,len(vars)):
            var = vars[j]
            dataset = np.zeros((grid_size,grid_size,grid_size))
            for i in range(0,num_patches):
                data_file = subdirectory + '/' + "processor_cluster.0000" + str(i) + ".samrai"
                data_h5file = h5py.File(data_file,'r')
                lower = patch_extents[i]["lower"]
                upper = patch_extents[i]["upper"]
                x_lower,x_upper = int(lower[0]),int(upper[0])+1
                y_lower,y_upper = int(lower[1]),int(upper[1])+1
                z_lower,z_upper = int(lower[2]),int(upper[2])+1
                dataset_name = "/processor.0000" + str(i) + "/level.00000" + "/patch.0000" + str(i) + "/" + var
                #! think about the potential reasons to mark order="F"
                # using numpy slicing to accelerate relative to for loop
                dataset[x_lower:x_upper,y_lower:y_upper,z_lower:z_upper] = data_h5file[dataset_name][...].reshape((x_upper-x_lower,y_upper-y_lower,z_upper-z_lower),order='F')
            data[var] = dataset
        return data
    out_step = treat_string(str(time*step_interval))
    subdirectory = main + "visit_dump." + out_step
    return extract_data(subdirectory,size)       

directory,step,start,threshold = sys.argv[1],int(sys.argv[2]),int(sys.argv[3]), int(sys.argv[4])
for i in range(start,threshold):
    data_all = manipulate_data(lattice_dim,i,step,directory)
    data_all = sum_modes(XYZ(),lattice_dim,(gw_fft(data_all)))
    density_spectrum = gw_spectrum(data_all,lattice_dim)
    write_file = 'OutData/' + str(i*step) + ".txt"
    f = open(write_file,'w')
    f.write(str(density_spectrum))
    f.close()
