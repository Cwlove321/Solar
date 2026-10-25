#!/usr/bin/env python3
"""Build an isolated, validated event case from preprocessed magnetogram/axis."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import shutil
import tempfile
from pathlib import Path

from atmosphere_config import FORTRAN_KEYS, normalize_atmosphere
from case_inputs import inspect_inputs, integer, number
from configure_case import replace_key_text
from refresh_par_readmes import refresh_all

ROOT = Path(__file__).resolve().parents[1]
STAGES = ("01_pfss","01_pfss_global_ct","02_relax","03_insert","04_fullmhd")
PAR_PATHS = [Path(s)/"amrvac.par" for s in STAGES] + [Path("04_fullmhd/profiles")/n for n in ("amrvac_full_thermo.par","amrvac_norad_with_heating.par")]
DYNAMIC_ITFIXGRID = 2147483647  # fixgrid uses it>=itfixgrid; negative values freeze immediately.


def _bool(value, name):
    if not isinstance(value,bool):
        raise ValueError(f"{name} must be a JSON boolean")
    return value


def _json(path: Path, value) -> None:
    path.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")


def _science(config: dict) -> dict:
    if not isinstance(config,dict):
        raise ValueError("settings must be a JSON object")
    allowed = {"schema_version","boundary_dir","path_dat","atmosphere","grid","rbsl","filament","insert_stage","stage_times","pfss_stage","fullmhd_source","thermal_profile","amr_mode","global_pfss","allow_reference_time_mismatch"}
    if set(config)-allowed:
        raise ValueError(f"unknown settings keys: {sorted(set(config)-allowed)}")
    resolved = dict(config)
    resolved["allow_reference_time_mismatch"] = _bool(config.get("allow_reference_time_mismatch",False),"allow_reference_time_mismatch")
    resolved["atmosphere"] = normalize_atmosphere(config["atmosphere"]) if "atmosphere" in config else normalize_atmosphere()
    for name, options, default in (("pfss_stage",("01_pfss","01_pfss_global_ct"),"01_pfss"),("fullmhd_source",("insert","relax"),"insert"),("thermal_profile",("default","full_thermo","norad_with_heating"),"default"),("amr_mode",("fixed","dynamic"),"fixed")):
        resolved[name] = config.get(name,default)
        if resolved[name] not in options:
            raise ValueError(f"{name} must be one of {options}")
    rbsl = dict(config["rbsl"])
    if set(rbsl)-{"bflux_mx","radius_mm","positive_helicity","flip_axis","npoints"}:
        raise ValueError("unknown rbsl settings")
    for key in ("bflux_mx","radius_mm"):
        rbsl[key] = number(rbsl[key],f"rbsl.{key}",positive=True)
    for key in ("positive_helicity","flip_axis"):
        rbsl[key] = _bool(rbsl[key],f"rbsl.{key}")
    resolved["rbsl"] = rbsl
    filament = dict(config["filament"])
    if set(filament)-{"method","density_factor","core_fraction","minimum_height_mm","liakh_mass_loading"}:
        raise ValueError("unknown filament settings")
    method = filament.get("method","none")
    if method not in ("none","instant","liakh_source"):
        raise ValueError("filament.method must be none, instant or liakh_source")
    filament["method"] = method
    for key,default in (("density_factor",20.),("core_fraction",.6),("minimum_height_mm",5.)):
        filament[key] = number(filament.get(key,default),f"filament.{key}",positive=True)
    if filament["density_factor"] < 1 or not 0 < filament["core_fraction"] < 1:
        raise ValueError("filament density_factor must be >=1 and core_fraction must be in (0,1)")
    liakh = dict(filament.get("liakh_mass_loading",{}))
    if set(liakh)-{"density_contrast","start_time_s","load_time_s","sigma_horizontal_mm","sigma_radial_mm","minimum_height_mm"}:
        raise ValueError("unknown Liakh settings")
    for key,default in (("density_contrast",30.),("start_time_s",0.),("load_time_s",100.),("sigma_horizontal_mm",3.6),("sigma_radial_mm",4.8),("minimum_height_mm",5.)):
        liakh[key] = number(liakh.get(key,default),f"filament.liakh_mass_loading.{key}",positive=key!="start_time_s")
    if liakh["start_time_s"] < 0:
        raise ValueError("Liakh start_time_s must be >=0")
    filament["liakh_mass_loading"] = liakh
    if resolved["fullmhd_source"] == "relax" and method != "none":
        raise ValueError("fullmhd_source=relax requires filament.method=none (no inserted rope)")
    resolved["filament"] = filament
    times = {}
    defaults = {"02_relax":{"time_max_code":1.,"save_interval_code":.05,"reset_time":True,"autoconvert":False},"03_insert":{"time_max_code":0.,"save_interval_code":.01,"reset_time":True,"autoconvert":False},"04_fullmhd":{"time_max_code":.3,"save_interval_code":.005,"reset_time":True,"autoconvert":True}}
    overrides = config.get("stage_times",{})
    if set(overrides)-set(defaults):
        raise ValueError("stage_times may contain only 02_relax, 03_insert and 04_fullmhd")
    for stage,values in defaults.items():
        entry = {**values,**(config.get("insert_stage",{}) if stage=="03_insert" else {}),**overrides.get(stage,{})}
        if set(entry)-set(values):
            raise ValueError(f"unknown stage time settings for {stage}")
        entry["time_max_code"] = number(entry["time_max_code"],f"{stage}.time_max_code")
        entry["save_interval_code"] = number(entry["save_interval_code"],f"{stage}.save_interval_code",positive=True)
        if entry["time_max_code"] < 0 or (stage!="03_insert" and entry["time_max_code"]==0):
            raise ValueError(f"{stage}.time_max_code must be positive (zero is allowed only for insertion)")
        for key in ("reset_time","autoconvert"):
            entry[key] = _bool(entry[key],f"{stage}.{key}")
        times[stage] = entry
    resolved["stage_times"] = times
    resolved["stages"] = [resolved["pfss_stage"],"02_relax"]+(["03_insert"] if resolved["fullmhd_source"]=="insert" else [])+["04_fullmhd"]
    previous_end = 0.
    for stage in resolved["stages"][1:]:
        entry = times[stage]
        start = 0. if entry["reset_time"] else previous_end
        target = entry["time_max_code"]
        if target < start-1e-12:
            raise ValueError(f"{stage}.time_max_code={target} is below its restart start time {start}; increase the absolute target or set reset_time=true")
        if stage=="03_insert" and method=="liakh_source":
            start_s,target_s = start*5972.477,target*5972.477
            if liakh["start_time_s"] < start_s-1e-4 or liakh["start_time_s"]+liakh["load_time_s"] > target_s+1e-4:
                raise ValueError(f"03_insert Liakh loading window must fit completely inside the actual stage time [{start_s:.9g},{target_s:.9g}] s; reset_time=false preserves the restart clock (unit_time approximately 5972.477 s)")
        previous_end = target
    return resolved


def _copy_template(target: Path) -> None:
    ignored = shutil.ignore_patterns(".DS_Store","__pycache__","*.pyc",".git","data","log","logs","amrvac","restart.dat","path.dat","Bphot","potential_boundary","makefile","mod_usr.f","*.o","*.mod")
    for stage in STAGES:
        shutil.copytree(ROOT/stage,target/stage,ignore=ignored)
    for directory in ("scripts","Preprocess"):
        shutil.copytree(ROOT/directory,target/directory,ignore=ignored)
    for name in ("README.md","QUICKSTART.md","DEPENDENCIES.md","VALIDATION.md","workflow.py","case.example.json","case.delivery.example.json","machine.example.json","machine.local.example.json","machine.lsf.example.json","tested-environment.txt",".gitignore"):
        if (ROOT/name).is_file():
            shutil.copy2(ROOT/name,target/name)
    for requirements in ROOT.glob("requirements*.txt"):
        shutil.copy2(requirements,target/requirements.name)
    _copy_locked_sources(target)


def _copy_locked_sources(target: Path) -> None:
    """Copy the source manifest only, never solver build/cache directories.

    The AMRVAC root is deliberately not passed through _copy_template's
    ignore pattern for the stage executable named ``amrvac``.
    """
    lock_path = ROOT/"vendor/dependencies.lock.json"
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    (target/"vendor").mkdir()
    shutil.copy2(lock_path,target/"vendor/dependencies.lock.json")
    shutil.copy2(ROOT/"vendor/SOURCE_PROVENANCE.md",target/"vendor/SOURCE_PROVENANCE.md")
    for package in lock["packages"].values():
        package_path = Path(package["path"])
        if package_path.is_absolute() or ".." in package_path.parts or package_path.parts[0] != "vendor":
            raise ValueError("dependency lock package must point inside vendor/")
        for relative,expected_hash in package["file_sha256"].items():
            relative_path = Path(relative)
            if relative_path.is_absolute() or ".." in relative_path.parts or "__pycache__" in relative_path.parts or relative_path.suffix in {".pyc",".o",".mod",".a",".so",".dylib"}:
                raise ValueError(f"dependency lock must list source files only: {relative}")
            source = ROOT/package_path/relative_path
            if source.is_symlink() or not source.is_file():
                raise ValueError(f"locked dependency source is missing/nonregular: {source}")
            if hashlib.sha256(source.read_bytes()).hexdigest() != expected_hash:
                raise ValueError(f"locked dependency source changed: {source}; verify the source snapshot before creating a case")
            destination = target/package_path/relative_path
            destination.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(source,destination)


def _configure(target: Path, resolved: dict) -> None:
    grid = resolved["grid"]
    rope = resolved["rbsl"]
    loaded = resolved["filament"]
    for relative in PAR_PATHS:
        path = target/relative
        text = path.read_text(encoding="utf-8")
        def put(key,value):
            nonlocal text
            text = replace_key_text(text,path,key,value)
        for key,ft in FORTRAN_KEYS.items():
            put(ft,resolved["atmosphere"][key])
        for key in ("domain_nx1","domain_nx2","domain_nx3","xprobmin1","xprobmax1","xprobmin2","xprobmax2","xprobmin3","xprobmax3","refine_max_level"):
            put(key,grid[key])
        put("qstretch_baselevel(1)",grid["qstretch_baselevel1"])
        put("itfixgrid",0 if resolved["amr_mode"]=="fixed" else DYNAMIC_ITFIXGRID)
        for key,value in {"Bflux":rope["bflux_mx"],"Rmfr":rope["radius_mm"],"I_Helix":rope["positive_helicity"],"npoints":rope["npoints"],"flip_axis":rope["flip_axis"]}.items():
            put(key,value)
        stage = relative.parts[0]
        put("do_rbsl",resolved["fullmhd_source"]=="insert" and stage in ("01_pfss","01_pfss_global_ct","03_insert"))
        if stage == "01_pfss_global_ct":
            for key,value in resolved["global_pfss"].items():
                if key in ("patch_x0","patch_y0","expected_nphi","expected_ntheta"):
                    put("global_"+key,value)
            put("global_pfss_lmax",resolved["global_pfss"]["lmax"])
            put("global_pfss_r0",resolved["global_pfss"]["r0"])
            put("global_pfss_rs",resolved["global_pfss"]["rs"])
        if stage in resolved["stage_times"]:
            timing = resolved["stage_times"][stage]
            for key,value in {"time_max":timing["time_max_code"],"dtsave_dat":timing["save_interval_code"],"reset_time":timing["reset_time"],"autoconvert":timing["autoconvert"]}.items():
                put(key,value)
        if stage=="03_insert":
            put("firstprocess",resolved["fullmhd_source"]=="insert")
            put("load_filament_density",loaded["method"]=="instant")
            for key,value in {"filament_density_factor":loaded["density_factor"],"filament_core_fraction":loaded["core_fraction"],"filament_min_height_Mm":loaded["minimum_height_mm"]}.items():
                put(key,value)
            liakh = loaded["liakh_mass_loading"]
            put("liakh_mass_loading",loaded["method"]=="liakh_source")
            for key,value in {"liakh_density_contrast":liakh["density_contrast"],"liakh_t_start_s":liakh["start_time_s"],"liakh_t_load_s":liakh["load_time_s"],"liakh_sigma_horizontal_Mm":liakh["sigma_horizontal_mm"],"liakh_sigma_radial_Mm":liakh["sigma_radial_mm"],"liakh_min_height_Mm":liakh["minimum_height_mm"]}.items():
                put(key,value)
        path.write_text(text,encoding="utf-8")
    profile = resolved["thermal_profile"]
    if profile != "default":
        shutil.copy2(target/"04_fullmhd/profiles"/f"amrvac_{profile}.par",target/"04_fullmhd/amrvac.par")


def create_case(boundary_dir: Path, path_dat: Path, settings_path: Path, output: Path) -> Path:
    """Create a new directory atomically; inputs/templates are never modified.

    Returns the resolved output Path. Existing output paths are refused. The
    resulting case has relative input paths and retains scientific settings,
    metadata and checksums in case.resolved.json and input_manifest.json.
    """
    boundary_dir,path_dat,settings_path,output = [Path(p).expanduser().resolve() for p in (boundary_dir,path_dat,settings_path,output)]
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    resolved = _science(json.loads(settings_path.read_text(encoding="utf-8")))
    output.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f".{output.name}.prepare-",dir=output.parent) as scratch:
        target = Path(scratch)/"case"
        target.mkdir()
        _copy_template(target)
        snapshot = target/"inputs/boundary"
        snapshot.mkdir(parents=True)
        for name in ("parameters.json","pfss_boundary.dat","Brll0000.fits","cutout.dat"):
            if (boundary_dir/name).is_file():
                shutil.copy2(boundary_dir/name,snapshot/name)
        axis_snapshot = target/"inputs/axis"
        axis_snapshot.mkdir()
        shutil.copy2(path_dat,axis_snapshot/"path.dat")
        for name in ("path_summary.json","path_lonlat_radius.txt"):
            if (path_dat.parent/name).is_file():
                shutil.copy2(path_dat.parent/name,axis_snapshot/name)
        shutil.copy2(settings_path,target/"inputs/settings.json")
        if not (snapshot/"Brll0000.fits").exists():
            spec = importlib.util.spec_from_file_location("rbsl_make_brll",ROOT/"Preprocess/make_pdfi_brll.py")
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            module.build_brll(snapshot,snapshot/"Brll0000.fits")
        grid_keys={"domain_nx1","domain_nx2","domain_nx3","xprobmin1","xprobmax1","xprobmin2","xprobmax2","xprobmin3","xprobmax3","refine_max_level","qstretch_baselevel1"}
        if set(resolved["grid"])-grid_keys:
            raise ValueError(f"unknown grid settings: {sorted(set(resolved['grid'])-grid_keys)}")
        checked = inspect_inputs(snapshot,axis_snapshot/"path.dat",resolved["grid"],resolved["rbsl"],allow_reference_time_mismatch=resolved["allow_reference_time_mismatch"])
        resolved["grid"] = checked["grid"]
        resolved["rbsl"]["npoints"] = checked["npoints"]
        gp = dict(resolved.get("global_pfss",{}))
        if set(gp)-{"lmax","r0","rs"}:
            raise ValueError("global_pfss settings may contain only lmax, r0 and rs; patch/shape are derived")
        gp.update(checked["global_pfss"])
        gp["lmax"] = integer(gp.get("lmax",360),"global_pfss.lmax")
        gp["r0"] = number(gp.get("r0",1.),"global_pfss.r0",positive=True)
        gp["rs"] = number(gp.get("rs",2.5),"global_pfss.rs",positive=True)
        if gp["rs"] <= max(gp["r0"],resolved["grid"]["xprobmax1"]):
            raise ValueError("global_pfss.rs must exceed both r0 and radial domain outer edge")
        if not math.isclose(gp["r0"],resolved["grid"]["xprobmin1"],abs_tol=1e-12):
            raise ValueError("global_pfss.r0 must match grid.xprobmin1")
        if gp["lmax"] > min(gp["expected_ntheta"]-1,gp["expected_nphi"]//2):
            raise ValueError("global_pfss.lmax exceeds global angular sampling; explicitly choose a compatible lmax in settings (the case builder never changes its scientific value automatically)")
        resolved["global_pfss"] = gp
        resolved["schema_version"] = 1
        resolved["boundary_dir"] = "inputs/boundary"
        resolved["path_dat"] = "inputs/axis/path.dat"
        resolved["input_origins"] = {"boundary":str(boundary_dir),"axis":str(path_dat),"settings":str(settings_path)}
        for stage in STAGES:
            destination = target/stage
            (destination/"Bphot").mkdir()
            (destination/"potential_boundary").mkdir()
            (destination/"data").mkdir()
            (destination/"log").mkdir()
            shutil.copy2(snapshot/"Brll0000.fits",destination/"Bphot/Brll0000.fits")
            shutil.copy2(snapshot/"pfss_boundary.dat",destination/"potential_boundary/pfss_boundary.dat")
            shutil.copy2(axis_snapshot/"path.dat",destination/"path.dat")
        _configure(target,resolved)
        refresh_all(target)
        manifest = {str(p.relative_to(target)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((target/"inputs").rglob("*")) if p.is_file()}
        _json(target/"input_manifest.json",{"algorithm":"sha256","files":manifest})
        resolved["input_manifest"] = "input_manifest.json"
        _json(target/"case.resolved.json",resolved)
        # Directory rename publishes a complete case only after every check.
        if output.exists():
            raise FileExistsError(f"output appeared during preparation: {output}")
        target.rename(output)
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--boundary",type=Path,required=True)
    parser.add_argument("--axis",type=Path,required=True)
    parser.add_argument("--settings",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args = parser.parse_args()
    try:
        print(create_case(args.boundary,args.axis,args.settings,args.output))
    except (ValueError,KeyError,OSError) as exc:
        parser.exit(2,f"case creation failed: {exc}\n")


if __name__=="__main__":
    main()
