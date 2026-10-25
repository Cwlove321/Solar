#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  PURPOSE:
# ;   Calculate the squashing factor Q at the photosphere or a cross section or a box volume, 
# ;   given a 3D magnetic field with Cartesian, uniform or stretched grids
# ; 
# ;   For details see:
# ;   Zhang, P., Chen, J.*, Liu, R. and Wang, C., 2022, FastQSL: A Fast Computation Method for Quasi-separatrix Layers. The Astrophysical Journal, 937, 26
# ;
# ;   FORTRAN procedures used
# ;   qfactor.f90
# ;   trace_bline.f90
# ;   trace_scott.f90
# ;
# ; ------COMPILATION 
# ;   For Linux and MacOS (either by ifort or gfortran): 
# ;      ifort -o qfactor.x qfactor.f90 -fopenmp -O3 -xHost -ipo
# ;   gfortran -o qfactor.x qfactor.f90 -fopenmp -Ofast -march=native
# ;
# ;   please specify the path of qfactor.x at the line of "spawn, 'qfactor.x' " in this file, 
# ;   or move qfactor.x to the $PATH (e.g. /usr/local/bin/) of the system
# ;
# ;   For Windows: (executing "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" in cmd first would be necessary)
# ;      ifort /o qfactor.exe qfactor.f90 /Qopenmp /O3 /QxHost /Qipo
# ;   please replace qfactor.x by qfactor.exe for the path specifing
# ;
# ;   -O3, -xHost, -ipo, -Ofast, -march=native are for a better efficiency;
# ;   -Ofast would be problematic for MacOS, then please substitutue -O3 for it; 
# ;
# ; INPUTS
# ;   Bx, By, Bz: 3D magnetic field, will be forcibly converted to float arrays while writing 'b3d.bin'
# ;   
# ; OPTIONAL INPUTS  
# ;
# ;   xa, ya, za:  coordinates of grids in 1D arrays. The size of xa, ya, za must be consistant with the size of the 3D magnetic field, 
# ;                and the values in x{yz}a should be increased by order
# ;   stretchFlag= keyword_set(xa) and keyword_set(ya) and keyword_set(za), which means bx, by, bz are on stretched grids, otherwise are on uniform grids; 
# ;                if invoked, the units of x{yz}reg, delta, length are same as x{yz}a, 
# ;                while the units of step and tol are the scale of a cell (a self-adaptive fashion that varying from cell to cell)
# ;
# ;   xreg, yreg, zreg: coordinates of the output region, in arrays of two elements, default is to include all available pixels at the photosphere.
# ;                 If stretchFlag=0B, the unit of coordinates is pixel, and the origin of the 3D field is [0, 0, 0].
# ;                 --- If (xreg[0] ne xreg[1]) AND (yreg[0] ne yreg[1]) AND (zreg[0] ne zreg[1]) AND NOT csFlag
# ;                     calculate Q in a box volume
# ;                     invoke vFlag
# ;                 --- If (xreg[0] eq xreg[1]) OR (yreg[0] eq yreg[1]) or (zreg[0] eq zreg[1])
# ;                     calculate Q in a cross section perpendicular to X or Y or Z axis
# ;                     invoke cFlag 
# ;                 --- If zreg=[0, 0], calculate Q at the photosphere
# ;                     invoke z0Flag 
# ;                 --- If csFlag is set, see below
# ;                     invoke cFlag
# ;   
# ;   csFlag:     to calculate Q in a cross section defined by three points; default is 0
# ;               point0=[xreg[0],yreg[0],zreg[0]] ; the origin of output
# ;               point1=[xreg[1],yreg[1],zreg[0]] ; point0 -> point1, first axis
# ;               point2=[xreg[0],yreg[0],zreg[1]] ; point0 -> point2, second axis
# ;               x{yz}reg[0] does not have to be smaller than x{yz}reg[1] in this case
# ;
# ;   factor:     to bloat up the original resolution, i.e. grid spacing of output = 1/factor; default is 4 
# ;
# ;   delta:      grid spacing of output; default is 1/factor (when stretchFlag=0B); 
# ;               if it is set, the keyword factor will be ignored
# ;
# ;   RK4Flag:    to trace bline by RK4; default is 0B (RKF45)
# ;   			   
# ;   step:       step size in tracing field lines for RK4; default=1.0
# ;
# ;   tol:        tolerance of a step in RKF45; default is 10.^(-4)
# ;
# ;   maxsteps:   maxium steps for stracing a field line at one direction; default is 4*(nx+ny+nz)/step; 
# ;               if highly twisted field lines exist and traced by RK4, this value should be larger; suggested by Jiang, Chaowei
# ;
# ;   scottFlag:  calculate Q and Q_perp by the method of Scott_2017_ApJ_848_117; default is 0B (method 3 of Pariat_2012_A&A_541_A78)
# ;		
# ;   twistFlag:  to calculate twist number Tw; see Liu_2016_ApJ_818_148; default is 0
# ;
# ;   curlB_out:  to save curlB at odir+'curlB.sav'; curlBx, curlBy, curlBz have same dimensions as Bx, By, Bz; default is 0
# ;
# ;   odir:       directory to save the results
# ;		
# ;   fstr:       filename of the results; file_sav=odir+fstr+'.sav'
# ;
# ;   nbridges:   number of processors to engage; default is !CPU.HW_NCPU - 2
# ; 
# ;   no_preview: don't produce PNG images for preview; default is 0B
# ;
# ;   RAMtmp:     use RAM (/dev/shm/tmp/) to speed up the data transmission, please only invoke it in Linux; default is 0B; 
# ;		if invoked, please run only one task of qfactor.pro simultaneously on one machine
# ;
# ;   tmpB:       apply temporary() to Bx, By, Bz to reduce the memory occupation of 3D magnetic field in IDL after writing b3d.bin to tmp_dir; default is 0B
# ;
# ;  memory occupation in qfactor.x: 
# ;       3D magnetic field + some 2D arrays + 3D curlB field (same as the occupation of the 3D magnetic field, when twistFlag=1B) +
# ;       grad_unit_vec_Bfield (3 times as the occupation of the 3D magnetic field)
# ;
# ;
# ; OUTPUTS:
# ;   q/qcs/q3d:  squashing factor 
# ;  
# ;   q_perp/q_perp3d:     q_perp in Titov 2007, ApJ, 660, 863
# ;
# ;   slogq/slogq_perp:    sign(Bz) x log_10(Q/Q_\perp); Calculated with all the field lines
# ;   slogq_orig/slogq_perp_orig:  only include field lines with both footpoints on the photoshere
# ;  
# ;   twist/twist3d:       twist number = \int \mu_0 J \cdot B /(4\pi) B^2 ds
# ;
# ;   rboundary/rsboundary/reboundary: (r:remote, s:start point, e:end point)
# ;                nature of the ends of field lines, see 'subroutine vp_rboundary' in trace_bline.f90
# ;                0 - inside
# ;                1 - end at zmin
# ;                2 - end at zmax
# ;                3 - end at ymin
# ;                4 - end at ymax
# ;                5 - end at xmin
# ;                6 - end at xmax
# ;                7 - others
# ;   rboundary3d=rsboundary3d+8*reboundary3d, for saving storage; So if rboundary3d[i, j, k] eq 9B, both two mapping surfaces of q3d[i, j, k] are the photosphere
# ;
# ;   rF/rsF/reF: coordinates of mapping points (F:foot); suggested by Jiang, Chaowei
# ; 
# ;   length: length of field lines
# ;
# ;   Bnr: abs(Bn_local/Bn_target) at the photosphere
import os
import numpy as np
import subprocess
import shutil
import time
import math
from pathlib import Path
from typing import Optional, Union

__version__ = "0.1.0"


def _find_qfactor_executable() -> Path:
    """Return the configured or bundled FastQSL CPU executable."""
    configured = os.environ.get("FASTQSL_EXECUTABLE")
    executable = Path(configured).expanduser() if configured else Path(__file__).parent / "bin" / "qfactor.x"
    executable = executable.resolve()
    if not executable.is_file():
        raise FileNotFoundError(
            f"FastQSL executable not found: {executable}. "
            "Run Packages/qsl/scripts/build_qfactor.sh or set FASTQSL_EXECUTABLE."
        )
    if not os.access(executable, os.X_OK):
        raise PermissionError(f"FastQSL executable is not executable: {executable}")
    return executable

def qfactor(
    Bx: np.ndarray,
    By: np.ndarray,
    Bz: np.ndarray,
    xa: Optional[np.ndarray] = None,
    ya: Optional[np.ndarray] = None,
    za: Optional[np.ndarray] = None,
    xreg: Optional[Union[list, np.ndarray]] = None,
    yreg: Optional[Union[list, np.ndarray]] = None,
    zreg: Optional[Union[list, np.ndarray]] = None,
    csFlag: bool = False,
    factor: int = 4,
    delta: Optional[float] = None,
    RK4Flag: bool = False,
    step: float = 1.0,
    tol: float = 1e-4,
    maxsteps: Optional[int] = None,
    scottFlag: bool = False,
    twistFlag: bool = False,
    save_curlB: bool = False,
    odir: Optional[str] = None,
    fstr: Optional[str] = None,
    nbridges: Optional[int] = None,
    no_preview: bool = False,
    tmpB: bool = False,
    RAMtmp: bool = False,
    compress: bool = False
):
    """
    Python 接口函数，用于调用 Fortran 程序 qfactor.x 计算 Q、Q_perp 等。
    参数含义与原 qfactor.pro 相同或相似；内部流程基本对应原 IDL 代码逻辑。
    """

    # -----------------------------------------------------------------------------
    # 0) 参数检查
    # -----------------------------------------------------------------------------
    # 形状检查：必须是三维数组 (nx, ny, nz)
    # 注意：如果您希望在 Fortran 中把第一维当作 x、第二维当作 y、第三维当作 z，
    #       而且 Fortran/IDL 是列优先，那么最好将 Python 中 Bx,By,Bz 先转置成 (nz,ny,nx)。
    #       下文我们会在写入前再做一次 transpose((2,1,0)) 来达到类似目的。

    
    if Bx.ndim != 3 or By.ndim != 3 or Bz.ndim != 3:
        raise ValueError("Bx, By, Bz 必须是 3 维数组 (nx, ny, nz)!")
    nz, ny, nx = Bx.shape
    if By.shape != (nz, ny, nx) or Bz.shape != (nz, ny, nx):
        raise ValueError("Bx, By, Bz 的维度必须完全一致！")

    # 检查是否有 (Bx, By, Bz) 全为 0 的点过多，或是否存在 NaN / Inf
    all_zero_mask = (Bx == 0.0) & (By == 0.0) & (Bz == 0.0)
    n_nulls = np.count_nonzero(all_zero_mask)
    if n_nulls > 3:
        raise ValueError(f"网格中存在 {n_nulls} 个以上的 (Bx,By,Bz) = (0,0,0) 点，可能导致场线追踪问题！")

    if not np.all(np.isfinite(Bx)) or not np.all(np.isfinite(By)) or not np.all(np.isfinite(Bz)):
        raise ValueError("网格中存在无穷大或 NaN 的 B 值，可能导致场线追踪问题！")

    # stretchFlag 是否开启
    stretchFlag = (xa is not None) and (ya is not None) and (za is not None)

    # xreg, yreg, zreg 默认值
    if xreg is None:
        if stretchFlag:
            xreg = [xa[0], xa[-1]]
        else:
            xreg = [0, nx - 1]
    if yreg is None:
        if stretchFlag:
            yreg = [ya[0], ya[-1]]
        else:
            yreg = [0, ny - 1]
    if zreg is None:
        if stretchFlag:
            zreg = [za[0], za[0]]  # 默认只取一个截面
        else:
            zreg = [0, 0]          # 默认 z=0

    # 转成 float
    xreg = np.array(xreg, dtype=np.float32)
    yreg = np.array(yreg, dtype=np.float32)
    zreg = np.array(zreg, dtype=np.float32)

    # csFlag 强制为 bool
    csFlag = bool(csFlag)

    # 线程数 nbridges
    if nbridges is None:
        import multiprocessing
        n_cpu = multiprocessing.cpu_count()
        nbridges = max(1, n_cpu - 2)

    # 如果没指定 delta，则默认 = 1/factor (当非 stretched 网格)
    if delta is None:
        if not stretchFlag:
            delta = 1.0 / factor
        else:
            delta = (xa[-1] - xa[0]) / ((nx - 1.0) * factor)
    delta = float(delta)

    # 若未指定 maxsteps，则默认 4*(nx+ny+nz)/step
    if maxsteps is None:
        maxsteps = int(4*(nx+ny+nz) / step)

    # -----------------------------------------------------------------------------
    # 1) 确定输出目录
    # -----------------------------------------------------------------------------
    cdir = os.getcwd()
    if odir is None:
        odir = os.path.join(cdir, "qfactor")
    odir = os.path.abspath(os.path.expanduser(os.fspath(odir)))
    if not os.path.exists(odir):
        os.makedirs(odir)

    # -----------------------------------------------------------------------------
    # 2) 准备临时目录，用于与 qfactor.x 通信
    # -----------------------------------------------------------------------------
    if RAMtmp:
        tmp_dir = "/dev/shm/FastQSLtmp/"
    else:
        tmp_dir = os.path.join(odir, "FastQSLtmp")

    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)
    os.makedirs(tmp_dir)

    # -----------------------------------------------------------------------------
    # 3) 写出 head.txt 等输入信息给 Fortran
    # -----------------------------------------------------------------------------
    head_txt_path = os.path.join(tmp_dir, "head.txt")
    with open(head_txt_path, "w") as f:
        # 顺序与 IDL 相同
        f.write(f"{nx} {ny} {nz} {nbridges} {delta} {maxsteps}\n")
        f.write(f"{xreg[0]} {xreg[1]} {yreg[0]} {yreg[1]} {zreg[0]} {zreg[1]} {step} {tol}\n")
        f.write(f"{int(twistFlag)} {int(RK4Flag)} {int(scottFlag)} {int(csFlag)} {int(save_curlB)}\n")

    # -----------------------------------------------------------------------------
    # 4) 写出 b3d.bin (Bx, By, Bz)
    #    由于 Fortran/IDL 列优先且往往把下标(1,1,2)当作( x=1, y=1, z=2 )，
    #    这里我们转置 Python 数组 (nx,ny,nz)->(nz,ny,nx)，再 asfortranarray() 写出
    # -----------------------------------------------------------------------------
    b3d_path = os.path.join(tmp_dir, "b3d.bin")
    with open(b3d_path, "wb") as f:

        np.asfortranarray(Bx, dtype=np.float32).tofile(f)
        np.asfortranarray(By, dtype=np.float32).tofile(f)
        np.asfortranarray(Bz, dtype=np.float32).tofile(f)

    # 若是 stretched 网格，也写 xa, ya, za (不做转置，因为它们是一维即可)
    if stretchFlag:
        np.asfortranarray(xa, dtype=np.float32).tofile(os.path.join(tmp_dir, "xa.bin"))
        np.asfortranarray(ya, dtype=np.float32).tofile(os.path.join(tmp_dir, "ya.bin"))
        np.asfortranarray(za, dtype=np.float32).tofile(os.path.join(tmp_dir, "za.bin"))

    # -----------------------------------------------------------------------------
    # 5) 调用外部 Fortran 程序 qfactor.x
    # -----------------------------------------------------------------------------
    t0 = time.time()
    executable = _find_qfactor_executable()
    ret = subprocess.run(
        [os.fspath(executable)],
        cwd=tmp_dir,
        capture_output=True,
        text=True,
        check=False,
    )
    t1 = time.time()
    if ret.returncode != 0:
        print("stdout:", ret.stdout)
        print("stderr:", ret.stderr)
        raise RuntimeError("调用 qfactor.x 失败，请检查是否编译并加入 PATH 或改为绝对路径！")

    print(f"qfactor.x 运行完成，耗时约 {t1 - t0:.2f} 秒")

    # -----------------------------------------------------------------------------
    # 6) 读取 tail.txt，判断实际输出类型
    # -----------------------------------------------------------------------------
    tail_path = os.path.join(tmp_dir, "tail.txt")
    with open(tail_path, "r") as f:
        lines = f.read().strip().split()
        qx = int(lines[0])
        qy = int(lines[1])
        qz = int(lines[2])
        q1 = int(lines[3])
        q2 = int(lines[4])
        z0Flag_ = int(lines[5])
        cFlag_  = int(lines[6])
        vFlag_  = int(lines[7])

    # -----------------------------------------------------------------------------
    # 7) 根据情况读取 Fortran 输出文件
    # -----------------------------------------------------------------------------
    results = {}
    if fstr is None:
        # 与 IDL 类似，根据情形确定前缀
        if z0Flag_ == 1:
            head_str = "qz0"
        elif cFlag_ == 1:
            head_str = "qcs"
        elif vFlag_ == 1:
            head_str = "q3d"
        else:
            head_str = "qfactor"
        delta_str = f"_delta{delta:.2f}"
        fstr = head_str + delta_str

    file_sav = os.path.join(odir, fstr + ".npz")

    # 如果要保存 curlB
    if save_curlB:
        curlB_path = os.path.join(tmp_dir, "curlB.bin")
        if os.path.isfile(curlB_path):
            curlB = np.fromfile(curlB_path, dtype=np.float32)
            # Fortran 里 curlB shape = (3, nx, ny, nz)
            # 读进后 reshape 一下
            curlB = curlB.reshape((3, nx, ny, nz), order='F')
            results["curlB"] = curlB
        np.savez_compressed(file_sav, **results)
        shutil.rmtree(tmp_dir)
        print(f"curlB 已保存到 {file_sav}")
        return results

    # z0Flag_ => photosphere
    if z0Flag_ == 1:
        path_q0 = os.path.join(tmp_dir, "qfactor0.bin")
        with open(path_q0, "rb") as f:
            count_q = q1 * q2
            q_arr = np.fromfile(f, dtype=np.float32, count=count_q).reshape((q1, q2), order='F')
            length_arr = np.fromfile(f, dtype=np.float32, count=count_q).reshape((q1, q2), order='F')
            bnr_arr = np.fromfile(f, dtype=np.float32, count=count_q).reshape((q1, q2), order='F')
            reF_arr = np.fromfile(f, dtype=np.float32, count=3*count_q).reshape((3, q1, q2), order='F')
            # Fortran writes reboundary as INTEGER(KIND=1), i.e. one byte per cell.
            rboundary_arr = np.fromfile(f, dtype=np.int8, count=count_q).reshape((q1, q2), order='F')

        results["q"] = q_arr
        results["length"] = length_arr
        results["Bnr"] = bnr_arr
        results["rF"] = reF_arr
        results["rboundary"] = rboundary_arr

        slogq_path = os.path.join(tmp_dir, "slogq.bin")
        slogq_arr = np.fromfile(slogq_path, dtype=np.float32)
        slogq_arr = slogq_arr.reshape((q1, q2), order='F')
        results["slogq"] = slogq_arr

        if scottFlag:
            q_perp_path = os.path.join(tmp_dir, "q_perp.bin")
            q_perp_arr = np.fromfile(q_perp_path, dtype=np.float32)
            q_perp_arr = q_perp_arr.reshape((q1, q2), order='F')
            results["q_perp"] = q_perp_arr

            slogq_perp_path = os.path.join(tmp_dir, "slogq_perp.bin")
            slogq_perp_arr = np.fromfile(slogq_perp_path, dtype=np.float32)
            slogq_perp_arr = slogq_perp_arr.reshape((q1, q2), order='F')
            results["slogq_perp"] = slogq_perp_arr

        if twistFlag:
            twist_path = os.path.join(tmp_dir, "twist.bin")
            twist_arr = np.fromfile(twist_path, dtype=np.float32)
            twist_arr = twist_arr.reshape((q1, q2), order='F')
            results["twist"] = twist_arr

        np.savez_compressed(file_sav, **results)

    # cFlag_ => cross section
    elif cFlag_ == 1:
        path_qcs = os.path.join(tmp_dir, "qcs.bin")
        with open(path_qcs, "rb") as f:
            count_q = q1 * q2
            qcs_arr = np.fromfile(f, dtype=np.float32, count=count_q).reshape((q1, q2), order='F')
            length_arr = np.fromfile(f, dtype=np.float32, count=count_q).reshape((q1, q2), order='F')

            rsF_arr = np.fromfile(f, dtype=np.float32, count=3*count_q).reshape((3, q1, q2), order='F')
            reF_arr = np.fromfile(f, dtype=np.float32, count=3*count_q).reshape((3, q1, q2), order='F')

            # rsboundary, reboundary 若需要则继续读:
            # rsboundary_arr = ...
            # reboundary_arr = ...

        results["qcs"] = qcs_arr
        results["length"] = length_arr
        results["rsF"] = rsF_arr
        results["reF"] = reF_arr

        if scottFlag:
            q_perp_path = os.path.join(tmp_dir, "q_perp.bin")
            q_perp_arr = np.fromfile(q_perp_path, dtype=np.float32)
            q_perp_arr = q_perp_arr.reshape((q1, q2), order='F')
            results["q_perp"] = q_perp_arr

        if twistFlag:
            twist_path = os.path.join(tmp_dir, "twist.bin")
            twist_arr = np.fromfile(twist_path, dtype=np.float32)
            twist_arr = twist_arr.reshape((q1, q2), order='F')
            results["twist"] = twist_arr

        np.savez_compressed(file_sav, **results)

    # vFlag_ => 3D box
    elif vFlag_ == 1:
        q3d_path = os.path.join(tmp_dir, "q3d.bin")
        q3d_arr = np.fromfile(q3d_path, dtype=np.float32)
        q3d_arr = q3d_arr.reshape((qx, qy, qz), order='F')
        results["q3d"] = q3d_arr

        rb3d_path = os.path.join(tmp_dir, "rboundary3d.bin")
        rb3d_arr = np.fromfile(rb3d_path, dtype=np.float32)
        rb3d_arr = rb3d_arr.reshape((qx, qy, qz), order='F')
        results["rboundary3d"] = rb3d_arr

        if scottFlag:
            q_perp3d_path = os.path.join(tmp_dir, "q_perp3d.bin")
            q_perp3d_arr = np.fromfile(q_perp3d_path, dtype=np.float32)
            q_perp3d_arr = q_perp3d_arr.reshape((qx, qy, qz), order='F')
            results["q_perp3d"] = q_perp3d_arr

        if twistFlag:
            twist3d_path = os.path.join(tmp_dir, "twist3d.bin")
            twist3d_arr = np.fromfile(twist3d_path, dtype=np.float32)
            twist3d_arr = twist3d_arr.reshape((qx, qy, qz), order='F')
            results["twist3d"] = twist3d_arr

        # 可能还会有 slogq.bin
        slogq_path = os.path.join(tmp_dir, "slogq.bin")
        if os.path.isfile(slogq_path):
            slogq_arr = np.fromfile(slogq_path, dtype=np.float32)
            slogq_arr = slogq_arr.reshape((q1, q2), order='F')
            results["slogq_z0"] = slogq_arr

        np.savez_compressed(file_sav, **results)

    else:
        # 既不 z0Flag 也不 cFlag, vFlag => 未产生任何输出？
        print("既不是 z0Flag 也不是 cFlag / vFlag，可能没有产生任何输出文件？")
        np.savez_compressed(file_sav, **results)

    # -----------------------------------------------------------------------------
    # 8) 清理临时目录，并返回
    # -----------------------------------------------------------------------------
    shutil.rmtree(tmp_dir)
    
    print(f"结果已保存至 {file_sav}")
    return results
