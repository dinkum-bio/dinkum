import csv

from . import vfg, vfn


def output_biotapestry_csv(output_filename):
    def fill_ten(ls):
        return list(ls) + [""] * (10 - len(ls))

    fp = open(output_filename, "w", newline="")
    w = csv.writer(fp)

    w.writerow(fill_ten(["# Model Commands"]))
    w.writerow(fill_ten(["# Command Type", "Model Name", "Parent Model"]))
    w.writerow(fill_ten(["model", "root"]))

    base_name = "double neg"
    w.writerow(fill_ten(["model", base_name, "root"]))

    for tissue in vfn._tissues:
        w.writerow(fill_ten(["model", tissue.name, base_name]))

    w.writerow(fill_ten(["# Region Commands"]))
    w.writerow(
        fill_ten(["# Command Type", "Model Name", "Region Name", "Region Abbreviation"])
    )

    for n, tissue in enumerate(vfn._tissues):
        w.writerow(fill_ten(["region", base_name, f"a{n}", f"a{n}"]))

    for n, tissue in enumerate(vfn._tissues):
        w.writerow(fill_ten(["region", tissue.name, f"a{n}", f"a{n}"]))

    w.writerow(fill_ten(["# Standard Interactions"]))
    w.writerow(
        fill_ten(
            [
                "# Command Type",
                "Model Name",
                "Source Type",
                "Source Name",
                "Target Type",
                "Target Name",
                "Sign",
                "Source Region Abbrev",
                "Target Region Abbrev",
            ]
        )
    )

    for ix in vfg._rules:
        model_name = base_name
        source_type = "gene"
        source_region = "a0"
        target_region = "a0"
        for target, source, sign in ix.btp_autonomous_links():
            source_name = source.name
            target_name = target.name
            assert sign in ["positive", "negative"]

            w.writerow(
                [
                    "general",
                    model_name,
                    source_type,
                    source_name,
                    "gene",
                    target_name,
                    sign,
                    source_region,
                    target_region,
                ]
            )

    w.writerow(fill_ten(["# Signals"]))
    w.writerow(fill_ten(["# Standalone nodes"]))
    w.writerow(fill_ten(["# Interactions for Submodels"]))
    w.writerow(
        [
            "# Command Type",
            "Model Name",
            "Source Type",
            "Source Name",
            "Target Type",
            "Target Name",
            "Sign",
            "Source Region Abbrev",
            "Target Region Abbrev",
        ]
    )
    fp.close()
