import numpy as np

try:
    from .nlfff_core.cartesian import (
        multigrid_reduce_components,
        preprocess_cartesian_field,
    )
    from .tmf_core.dave4vm import do_dave4vm
    from .tmf_core.sharp_boundary import sub_data
except ImportError:  # Direct import from the Preprocess directory.
    from nlfff_core.cartesian import (
        multigrid_reduce_components,
        preprocess_cartesian_field,
    )
    from tmf_core.dave4vm import do_dave4vm
    from tmf_core.sharp_boundary import sub_data
from glob import glob
from tqdm import tqdm
import matplotlib.pyplot as plt


def reduce_components(bx, by, bz, *, level=1, preprocess=False):
    """Use the shared Cartesian reducer and optional boundary preprocessing."""

    bx, by, bz = multigrid_reduce_components(bx, by, bz, level=level)
    if preprocess:
        bx, by, bz, _ = preprocess_cartesian_field(bx, by, bz)
    return bx, by, bz


##存的是n-1个中间时刻的v，n-1个后面时刻的B，第二张B为tmf时间起点和nlfff外推的底图
def write_driven_boundary(path,outpath,x1,x2,y1,y2,level=1,prep=0,driven_fast = 12,window_size = 23,thr = 100,nlfff_level=2,nlfff_prep = 1):
    
    t_simulation = []  #t-第2张磁图
    real_time = []
    Vx,Vy,Vz,Bx,By,Bz = [],[],[],[],[],[]

    path_Br = glob(path + '*.Br.fits')
    path_Br.sort()
    path_Br = path_Br
    path_Bp = glob(path + '*.Bp.fits')
    path_Bp.sort()
    path_Bp = path_Bp
    path_Bt = glob(path + '*.Bt.fits')
    path_Bt.sort()
    path_Bt = path_Bt

    nz = len(path_Br)
    dx = 366.4125*2**(level-1)
    dy,dz = dx,dx
    
    bx_start,by_start,bz_start,t0 = sub_data(path_Bp[0],path_Bt[0],path_Br[0],x1,x2,y1,y2)
    bx_start,by_start,bz_start = reduce_components(
        bx_start, by_start, bz_start, level=level, preprocess=False
    )
    
    for n in tqdm(range(1,nz)):

        bx_stop,by_stop,bz_stop,t1 = sub_data(path_Bp[n],path_Bt[n],path_Br[n],x1,x2,y1,y2)
        if n==1:   #start nlfff
            t_start = t1
            write_nlfff_boundary(outpath,bx_stop.copy(),by_stop.copy(),bz_stop.copy(),nlfff_level = nlfff_level,nlfff_prep =nlfff_prep,xc=0,yc=0)
        
        bx_prep,by_prep,bz_prep = reduce_components(
            bx_stop, by_stop, bz_stop, level=level, preprocess=bool(prep)
        )
        Bx.append(bx_prep)
        By.append(by_prep)
        Bz.append(bz_prep)
        
        
        bx_stop,by_stop,bz_stop = reduce_components(
            bx_stop, by_stop, bz_stop, level=level, preprocess=False
        )
        real_time.append(t1)
        dt = (t1-t0).total_seconds()/driven_fast
        t_s = (t1-t_start).total_seconds()/driven_fast
        t_simulation.append(t_s)

        #Calling do_dave4vm
        magvm, vel4vm = do_dave4vm(dt,bx_stop,bx_start,by_stop,
                                   by_start,bz_stop,bz_start,dx,
                                   dy, window_size)
        bx_c,by_c,bz_c = magvm['bx'],magvm['by'],magvm['bz']
        b_total = np.sqrt(bx_c**2+by_c**2+bz_c**2)

        vx_c,vy_c,vz_c = np.zeros_like(bx_c),np.zeros_like(by_c),np.zeros_like(bz_c)
        vx_c[window_size:-window_size,window_size:-window_size] = vel4vm['U0'][window_size:-window_size,window_size:-window_size]
        vy_c[window_size:-window_size,window_size:-window_size] = vel4vm['V0'][window_size:-window_size,window_size:-window_size]
        vz_c[window_size:-window_size,window_size:-window_size] = vel4vm['W0'][window_size:-window_size,window_size:-window_size]   
        vx_c[b_total<thr] = 0   #let vmap clean
        vy_c[b_total<thr] = 0
        vz_c[b_total<thr] = 0
        vector_plot(bz_c,vx_c/12,vy_c/12,intervel=6,scale=15,time = t1, save_name= outpath+'fig/v_fig/'+('%04d'% n)+'.png')   #后一张的时间戳
        
        Vx.append(vx_c)
        Vy.append(vy_c)
        Vz.append(vz_c)
        
        bx_start,by_start,bz_start = bx_stop.copy(),by_stop.copy(),bz_stop.copy()
        t0 = t1

    ny,nx = bx_start.shape
    nz = len(Bz)
    
    # Write magnetic field data to binary file
    with open(outpath+'tmf/bd/B_boundary.dat', 'wb') as f:
        np.array([nx, ny, nz], dtype=np.int32).tofile(f)
        np.array([dx, dy, dz], dtype=np.float64).tofile(f)
        np.array([Bx, By, Bz], dtype=np.float64).tofile(f)
        np.array(t_simulation, dtype=np.float64).tofile(f)

    # Write velocity field data to binary file
    with open(outpath+'tmf/bd/V_boundary.dat', 'wb') as f:
        np.array([nx, ny, nz], dtype=np.int32).tofile(f)
        np.array([dx, dy, dz], dtype=np.float64).tofile(f)
        np.array([Vx, Vy, Vz], dtype=np.float64).tofile(f)
        np.array(t_simulation, dtype=np.float64).tofile(f)
    np.save(outpath+'tmf/bd/UTC.npy',np.array(real_time))

def write_nlfff_boundary(outpath,bx,by,bz,nlfff_level,nlfff_prep =0,xc=0,yc=0):

    arcsec2cm = 7.32825e7  # 1 arcsec in cm
    dx = 0.50*2**(nlfff_level-1)*arcsec2cm  # x spatial resolution of the preprocessed magnetic field
    dy = dx  # y spatial resolution
    bx,by,bz = reduce_components(
        bx, by, bz, level=nlfff_level, preprocess=bool(nlfff_prep)
    )
    ny, nx = bx.shape
        
    
    filename = outpath+'tmf/bd/potential_boundary.dat'
    with open(filename, 'wb') as f:
        np.array(nx, dtype=np.int32).tofile(f)
        np.array(ny, dtype=np.int32).tofile(f)
        np.array(xc, dtype=np.float64).tofile(f)
        np.array(yc, dtype=np.float64).tofile(f)
        np.array(dx, dtype=np.float64).tofile(f)
        np.array(dy, dtype=np.float64).tofile(f)
        np.array(bz, dtype=np.float64).tofile(f)
    
    bx_n,by_n,bz_n = np.zeros([ny+4,nx+4]),np.zeros([ny+4,nx+4]),np.zeros([ny+4,nx+4])
    bx_n[2:-2,2:-2] = bx
    by_n[2:-2,2:-2] = by
    bz_n[2:-2,2:-2] = bz
    
    filename = outpath+'tmf/bd/nlfff_boundary.dat'
    with open(filename, 'wb') as f:
        np.array(nx+4, dtype=np.int32).tofile(f)
        np.array(ny+4, dtype=np.int32).tofile(f)
        np.array(xc, dtype=np.float64).tofile(f)
        np.array(yc, dtype=np.float64).tofile(f)
        np.array(dx, dtype=np.float64).tofile(f)
        np.array(dy, dtype=np.float64).tofile(f)
        np.array(bx_n, dtype=np.float64).tofile(f)
        np.array(by_n, dtype=np.float64).tofile(f)
        np.array(bz_n, dtype=np.float64).tofile(f)

    
    x1 = xc - nx * dx / 2
    x2 = xc + nx* dx / 2
    y1 = yc - ny * dy / 2
    y2 = yc + ny * dy / 2

    # 打印范围信息
    print('x, y, and z range (10 Mm):')
    print('        xprobmin1=', f'{x1 * 1.e-9}d0')
    print('        xprobmax1=', f'{x2 * 1.e-9}d0')
    print('        xprobmin2=', f'{y1 * 1.e-9}d0')
    print('        xprobmax2=', f'{y2 * 1.e-9}d0')
    print('        xprobmin3=', f'{(0.0 * 1.e-9 + 0.1)}d0')  # to lift the domain 1 Mm above 0
    print('        xprobmax3=', f'{(y2 - y1) * 1.e-9}d0')

def vector_plot(bz,vx,vy,intervel,scale,time,save_name):
    mask = np.sqrt(vx**2+vy**2)*np.sign(bz)
    fig,ax = plt.subplots(figsize=(10,5))
    fig.dpi = 400
    ny,nx = np.shape(bz)
    x = np.arange(nx)
    y = np.arange(ny)
    X,Y = np.meshgrid(x,y)
    cax = ax.imshow(bz,cmap = 'gray',origin = 'lower',vmax = 1000,vmin = -1000)
    
    n = intervel #稀疏间隔
    
    vx_p = vx.copy()
    vx_p[mask<0.1] = np.nan
    vy_p = vy.copy()
    vy_p [mask<0.1] = np.nan
    
    quiver = ax.quiver(X[::n,::n],Y[::n,::n],vx_p[::n,::n],vy_p[::n,::n],scale = scale,
                      color = 'skyblue',width = 0.002,angles='xy',alpha = 1)
    
    vx_n = vx.copy()
    vx_n[mask>-0.1] = np.nan
    vy_n = vy.copy()
    vy_n [mask>-0.1] = np.nan
    
    quiver = ax.quiver(X[::n,::n],Y[::n,::n],vx_n[::n,::n],vy_n[::n,::n],scale =scale,
                      color = 'tab:pink',width = 0.002,angles='xy',alpha = 1)
    
    ax.quiverkey(quiver, 0.84, 0.95, 0.3, '{} km/s'.format(0.3),
                                coordinates='axes', labelpos='E', color='black') 
    ax.set_title(time,size=10)
    ax.set_axis_off() 
    plt.savefig(save_name)
    plt.close(fig)
