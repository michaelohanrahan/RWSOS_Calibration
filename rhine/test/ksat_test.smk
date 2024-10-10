from pathlib import Path

cwd = "/p/11209265-grade2023/wflow/RWSOS_Calibration/rhine/data/2-interim"
base = "/p/11209265-grade2023/wflow/RWSOS_Calibration/rhine"
mods = ["base", "brt", "rf"]
configfile: '/u/ohanrah/documents/RWSoS/RWSOS_Calibration/rhine/test/config.yaml'
rule all:
    input:
        expand(Path(cwd, "ksat_test_"+"{mod}", "output_scalar.nc"), mod=mods),
        Path(base, "data/5-visualization", "ksat_map_test", "interactive", "Peaktiming_Lobith_9_20150801-20180130.html")

rule wflow:
    input:
        cfg=Path(cwd, "ksat_test_"+"{mod}", "wflow_sbm.toml")
    params:
        project=Path(base, "bin")
    output:
        Path(cwd, "ksat_test_"+"{mod}", "output_scalar.nc")
    localrule: False
    group: 'wflow'
    threads: 1
    resources:
        mem_mb=8000
    shell:
        f"""julia --project="{{params.project}}" -t {{threads}} -e\
        "using Pkg;\
        Pkg.instantiate();\
        using Wflow;
        Wflow.run()"\
        {{input.cfg}}"""

rule done:
    input:
        Path(cwd, "ksat_test_"+"{mod}", "output_scalar.nc")
    output:
        done = Path(cwd, "ksat_test_"+"{mod}", "run.done")
    localrule: True
    shell:'''touch {output.done}'''
    

rule plot_comparison:
    input:
        files = expand(Path(cwd, "ksat_test_"+"{mod}", "output_scalar.nc"), mod=mods),
        dones = expand(Path(cwd, "ksat_test_"+"{mod}", "run.done"), mod=mods)
    params:
        script = Path(base, "src/post/plot_final_model.py"),
        workdir = Path(base),
        gaugetoplot = Path(cwd, "wflow_id_rhine.csv"),
        obs = Path(base, "data/1-external", "discharge_obs_hr_FORMAT_allvars_wflowid_0_to_727.nc"),
        starttime = "2015-08-01",
        endtime = "2018-01-30",
        outdir = Path(base, "data/5-visualization", "ksat_map_test")
    output:
        Path(base, "data/5-visualization", "ksat_map_test", "interactive", "Peaktiming_Lobith_9_20150801-20180130.html")
    localrule: True
    run:
        '''{params.script}'''