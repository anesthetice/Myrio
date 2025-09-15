import marimo

__generated_with = "0.15.2"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import polars as pl
    return mo, pl


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Bold Database Creation/Conversion

    1. "Raw" data downloaded as `.tsv` from `https://portal.boldsystems.org/result?query=Plantae[tax]` (31.08.2025)
    2. File opened and edited in `Libreoffice Calc`
        * Copy the following columns to a new tsv file (encode as a utf-8 for the love of god):
            * A: process_id
            * O: phylum
            * P: class
            * Q: order
            * R: family
            * U: genus
            * V: species
            * BO: nuc (sequence)
            * BS: marker_coode


    Data downloaded from
    """
    )
    return


@app.cell
def _(pl):
    raw_df = pl.read_csv(
        "./data/BOLD_Plantae_20250831_edited.tsv", separator="\t", truncate_ragged_lines=True
    )
    raw_df
    return (raw_df,)


@app.cell
def _(raw_df):
    for idx, (item, item_nb) in enumerate(
        zip(raw_df["marker_code"].unique(maintain_order=True), raw_df["marker_code"].unique_counts())
    ):
        print(f"{item} ({item_nb})", end=", ")
        if idx % 9 == 8:
            print("\n")
    return


@app.cell
def _(pl, raw_df):
    df = raw_df.filter(
        pl.col("nuc").str.len_bytes() > 20,
        pl.col("family").is_not_null(),
    )
    df
    return (df,)


@app.cell
def _(df, pl):
    df_ITS = df.filter(pl.col("marker_code").is_in(["ITS", "ITS1", "ITS2"]))
    df_matK = df.filter(pl.col("marker_code").is_in(["matK"]))
    df_rbcL = df.filter(pl.col("marker_code").is_in(["rbcL"]))
    df_trnH_psbA = df.filter(pl.col("marker_code").is_in(["trnH-psbA"]))
    return df_ITS, df_matK, df_rbcL, df_trnH_psbA


@app.cell
def _(df_ITS, df_matK, df_rbcL, df_trnH_psbA, pl):
    def process_and_save_df_to_fasta(df: pl.DataFrame, filepath: str):
        data = []
        for row in df.rows():
            pid, p, c, o, f, g, s, seq, _ = row
            # seq_lb = '\n'.join(seq[i:i+80] for i in range(0, len(seq), 80))
            string = f">BOLD_PROCESS_ID={pid}|tax={{"
            for clade, name in zip(["p", "c", "o", "f", "g", "s"], [p, c, o, f, g, s]):
                if name is not None:
                    if clade != "s":
                        string += f"{clade}:{name}, "
                    else:
                        string += f"{clade}:{name}"
            string += "}\n"
            string += seq.upper()
            data.append(string)

        data = "\n".join(data)
        file = open(filepath, "w", encoding="utf-8")
        file.write(data)
        file.close()


    process_and_save_df_to_fasta(df_ITS, "./data/BOLD_Plantae_20250831_ITS.fasta")
    process_and_save_df_to_fasta(df_matK, "./data/BOLD_Plantae_20250831_matK.fasta")
    process_and_save_df_to_fasta(df_rbcL, "./data/BOLD_Plantae_20250831_rbcL.fasta")
    process_and_save_df_to_fasta(df_trnH_psbA, "./data/BOLD_Plantae_20250831_trnH-psbA.fasta")
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
