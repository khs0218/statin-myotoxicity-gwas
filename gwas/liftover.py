def add_chr_prefix(x):
    x = str(x)
    if "chr" in x:
        return x
    else:
        return "chr" + x


def convert_and_update(row):
    chrom = row["CHR"]
    pos = row["POS"]
    result = lo.convert_coordinate(chrom, pos)

    if result:
        best_result = max(result, key=lambda x: x[3])
        row["Converted_Chromosome"] = best_result[0]
        row["Converted_Position"] = best_result[1]
        row["Converted_Strand"] = best_result[2]
        row["Conversion_Score"] = best_result[3]
    else:
        row["Not_Converted"] = "yes"

    return row


if __name__ == "__main__":
    import argparse
    import os
    import pandas as pd
    from pyliftover import LiftOver

    parser = argparse.ArgumentParser(description="Allofus Liftover")
    parser.add_argument("--data", type=str, default="")
    args = parser.parse_args()
    lo = LiftOver("hg38", "hg19")

    data_dir = "/result/allofus/"
    data_path = os.path.join(data_dir, args.data + ".txt")
    aou = pd.read_csv(data_path, sep="\t")
    aou["CHR"] = aou["CHR"].apply(add_chr_prefix)
    aou["Converted_Chromosome"] = pd.NA
    aou["Converted_Position"] = pd.NA
    aou["Converted_Strand"] = pd.NA
    aou["Conversion_Score"] = pd.NA
    aou["Not_Converted"] = pd.NA
    aou = aou.apply(convert_and_update, axis=1)
    count = aou["Not_Converted"].notna().sum()
    print(count)
    aou_filtered = aou[~aou["Converted_Chromosome"].isna()]
    aou_filtered = aou_filtered[
        aou_filtered["Converted_Chromosome"].isin([f"chr{i}" for i in range(1, 23)])
    ]
    aou_filtered["Converted_Chromosome"] = (
        aou_filtered["Converted_Chromosome"]
        .str.replace("chr", "", regex=False)
        .astype(int)
    )
    aou_filtered["Converted_Position"] = aou_filtered["Converted_Position"].astype(int)
    locus = aou_filtered[["Converted_Chromosome", "Converted_Position", "Allele1", "Allele2", "BETA", "VAR", "p.value"]]
    locus.to_csv(
        os.path.join(data_dir, "lifted", args.data + "_37.txt"), sep="\t", index=False
    )
