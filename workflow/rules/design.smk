"""
filename: design.smk

Description:
    Organize all fastq files within `config["fq_dir"]` according to 
    the following keywords:
    - wt (list): wt-1 wt-2 wt-3 ...
    - mut (list): mut-1 mut-2 mut-3 ...
    - fq_dir (str): /path/to/fq_files

    The fastq files within `config["fq_dir"]` are expected named according 
    to the following format:
    {sample_name}_{unit}_{read12}{ext}
    - only letters, digits, underscore, hyphen [A-Za-z0-9_.-]
    eg: RNAseq_mESC_DMSO_2h_rep1_1.fq.gz

Returns:
    Two files for further analysis:
    - results/design/samples.tsv
    - results/design/units.tsv

    # file: samples.tsv
    sample_name     group   treatment_1     treatment_2     jointly_handled
    RNAseq_mESC_DMSO_2h 1   wt              untreated       1
    RNAseq_mESC_dTAG_2h 1   mut             untreated       1

    # file: units.tsv
    sample_name     unit_name       fq1     fq2
    RNAseq_mESC_DMSO_2h rep1   RNAseq_mESC_DMSO_2h_rep1_1.fq.gz  RNAseq_mESC_DMSO_2h_rep1_2.fq.gz
    RNAseq_mESC_DMSO_2h rep2   RNAseq_mESC_DMSO_2h_rep2_1.fq.gz  RNAseq_mESC_DMSO_2h_rep2_2.fq.gz
    RNAseq_mESC_dTAG_2h rep1   RNAseq_mESC_dTAG_2h_rep1_1.fq.gz  RNAseq_mESC_dTAG_2h_rep1_2.fq.gz
    RNAseq_mESC_dTAG_2h rep2   RNAseq_mESC_dTAG_2h_rep2_1.fq.gz  RNAseq_mESC_dTAG_2h_rep2_2.fq.gz
"""

# import os
import sys
import warnings
from pathlib import Path

# load functions within scripts
scripts_dir = str(Path(workflow.basedir) / "scripts")
sys.path.append(scripts_dir)
from utils import fx_name, fx_name_parts, list_fq_dir, save_to_yaml


# Generate samples.tsv file
rule build_samples:
    output:
        config["samples"],
    run:
        # check wt and mut
        wt = config.get("wt", [])
        mut = config.get("mut", [])
        if len(wt) != len(mut) or len(wt) < 0:
            raise ValueError("design.wt and design.mut not valid")
            # all fastq files
        fq_lists = list_fq_dir(config["fq_dir"], paired_only=True)
        fq_names = [
            fx_name(i[0], fix_rep=True, fix_pe=True) for i in fq_lists
        ]  # all names
        fq_names = list(set(fq_names))  # all available names
        lines = []
        # for w in wt:
        for ia, w in enumerate(wt):
            lines.extend(
                [
                    "\t".join([i, f"{ia+1}", "wt", "untreated", "1"])
                    for i in fq_names
                    if w in i
                ]
            )
        for ib, m in enumerate(mut):
            lines.extend(
                [
                    "\t".join([i, f"{ib+1}", "mut", "untreated", "1"])
                    for i in fq_names
                    if m in i
                ]
            )
            # save design
        header = [
            "sample_name",
            "group",
            "treatment_1",
            "treatment_2",
            "jointly_handled",
        ]
        try:
            with open(output[0], "wt") as w:
                w.write("\t".join(header) + "\n")
                w.write("\n".join(lines) + "\n")
                # print("\n".join(lines))
        except Exception as e:
            warnings.warn(f"failed to save samples.tsv, {e}")


# Generate units.tsv file
## units could be: replicates (technical, biological)
## full_name: {sample}_{unit}_1.fq.gz
rule build_units:
    input:
        config["samples"],
    output:
        config["units"],
    run:
        # load sample_name
        df = pd.read_csv(input[0], sep="\t", dtype={"sample_name": str})
        samples = df.sample_name.values  # as list
        # read sample_name
        fq_lists = list_fq_dir(config["fq_dir"], paired_only=True)
        # output
        lines = []
        for sample in samples:
            for fq1, fq2 in fq_lists:
                if sample not in Path(fq1).name:
                    continue  #
                name, unit, read12, ext = fx_name_parts(fq1)
                ss = [name, unit, fq1, fq2]
                ss = list(map(str, ss))  # convert to list
                lines.append("\t".join(ss))
                # save units
        header = ["sample_name", "unit_name", "fq1", "fq2"]
        try:
            with open(output[0], "wt") as w:
                w.write("\t".join(header) + "\n")
                w.write("\n".join(lines) + "\n")
                # print("\n".join(lines))
        except Exception as e:
            warnings.warn(f"failed to save units.tsv, {e}")


