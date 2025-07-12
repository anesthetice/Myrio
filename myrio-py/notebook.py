import marimo

__generated_with = "0.14.10"
app = marimo.App(width="medium")


@app.cell
def _():
    from Bio import SeqIO
    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import polars as pl
    import scipy as sp
    from scipy import stats
    import seaborn as sns

    sns.set_theme()
    plt.rcParams["figure.dpi"] = 300
    plt.rcParams["axes.titlesize"] = 14
    plt.rcParams["axes.labelsize"] = 12
    plt.rcParams["xtick.labelsize"] = 9
    plt.rcParams["ytick.labelsize"] = 9
    plt.rcParams["legend.fontsize"] = 10
    return SeqIO, mo, np, pl, plt


@app.cell
def _(pl):
    data_df = pl.read_csv("data.csv")
    data_df
    return (data_df,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Quality Score Histogram""")
    return


@app.cell
def _(SeqIO, data_df, np, pl, plt):
    def _():
        df = data_df.filter(pl.col("filepath").str.ends_with(".fastq"))
        filepath = df.row(0)[6]

        counts = np.zeros(128)  # Phred quality score goes from 0 to 40
        for filepath in df["filepath"]:
            for record in SeqIO.parse(
                filepath, "fastq"
            ):  # fastq is equivalent to fastq-sanger → Phred quality score
                counts += np.bincount(record.letter_annotations["phred_quality"], minlength=128)
        freqs = counts / np.sum(counts)
    
        fig, ax1 = plt.subplots(1, 1, figsize=(9, 5))
        ax: plt.Axes = ax1
        ax.set_xlabel("Q-score")
        ax.set_ylabel("frequency")
        ax.set_xlim(0, 70)
    
        ax.bar(np.arange(128), freqs)
    

        return fig

    _()
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
