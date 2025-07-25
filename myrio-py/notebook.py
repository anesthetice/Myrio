import marimo

__generated_with = "0.14.12"
app = marimo.App(width="medium")


@app.cell
def _():
    from Bio import SeqIO
    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import polars as pl
    import pymsaviz
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
    return SeqIO, mo, np, pl, plt, pymsaviz


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

        counts = np.zeros(128)  # Phred quality score goes from 0 to 40
        for filepath in df["filepath"]:
            file_counts = np.zeros(128)
            # fastq is equivalent to fastq-sanger → Phred quality score
            for record in SeqIO.parse(filepath, "fastq"):
                file_counts += np.bincount(record.letter_annotations["phred_quality"], minlength=128)
            counts += file_counts
            print(f"{filepath} : {np.sum(file_counts * np.arange(128))/np.sum(file_counts)}")
        freqs = counts / np.sum(counts)

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 5))
        ax: plt.Axes = ax1
        ax.set_xlabel("Q-score")
        ax.set_ylabel("frequency")
        ax.set_xlim(0, 70)

        ax.bar(np.arange(128), freqs)

        ax: plt.Axes = ax2
        ax.set_xlabel("Q-score clean")
        ax.set_ylabel("frequency")
        ax.set_xlim(0, 41)

        freqs_clean = freqs[0:41] / np.sum(freqs[0:41])
        ax.bar(np.arange(41), freqs_clean)

        return freqs, freqs_clean, fig


    Q_score_freqs, Q_score_freqs_clean, _fig = _()
    Q_score_freqs_cummul_clean = np.cumulative_sum(Q_score_freqs_clean)

    print(Q_score_freqs_clean)
    print(Q_score_freqs_cummul_clean)

    _fig
    return


@app.cell
def _(SeqIO, data_df, np, plt):
    def _():
        filepath = data_df.row(0)[6]
        records = list(SeqIO.parse(filepath, "fastq"))
        record = records[2]
        qual = record.letter_annotations["phred_quality"]

        fig, ax1 = plt.subplots(1, 1, figsize=(20, 5))
        ax: plt.Axes = ax1
        ax.set_xlabel("idx")
        ax.set_ylabel("Q-score")

        ax.plot(np.arange(len(qual)), qual, ls="--", color="grey", mfc="blue", mec="blue", marker=".")
        ax.set_ylim(0, 40)
        ax.set_xlim(0, 856)

        return fig


    _()
    return


@app.cell
def _(SeqIO, data_df, np, pl, plt):
    def _():
        df = data_df.filter(pl.col("filepath").str.ends_with(".fastq"))

        counts = np.zeros(128)  # Phred quality score goes from 0 to 40
        for filepath in df["filepath"]:
            for record in SeqIO.parse(filepath, "fastq"):
                qual = record.letter_annotations["phred_quality"]
                blocks = [1]
                block_idx = 0
                block_main_qual = 0
                for q in qual[1:]:
                    if abs(block_main_qual - q) < 4:
                        blocks[block_idx] += 1
                    else:
                        block_main_qual = q
                        blocks.append(1)
                        block_idx += 1
                counts += np.bincount(blocks, minlength=128)[0:128]

        fig, ax1 = plt.subplots(1, 1, figsize=(9, 4))
        ax: plt.Axes = ax1
        ax.set_xlabel("frequency")
        ax.set_ylabel("block size")
        plt.gca().invert_yaxis()
        ax.barh(np.arange(30), (counts / np.sum(counts))[:30])

        freqs_safe = counts[:15] / np.sum(counts[:15])
        print(freqs_safe)
        freqs_safe_cummul = np.cumulative_sum(freqs_safe)
        print(freqs_safe_cummul)

        return fig


    _()
    return


@app.cell
def _():
    a = [1, 2, 3]  # todo: read from file
    return (a,)


@app.cell
def _(a, np, plt):
    def _():
        fig, ax1 = plt.subplots(1, 1, figsize=(9, 5))
        ax: plt.Axes = ax1
        ax.set_xlabel("Q-score")
        ax.set_ylabel("frequency")
        ax.set_xlim(0, 70)

        ax.bar(np.arange(41), np.bincount(a, minlength=41))

        fig.tight_layout()
        return fig


    _()
    return


@app.cell(hide_code=True)
def _():
    # fmt: off
    b = [27, 27, 29, 25, 27, 10, 11, 11, 7, 8, 10, 10, 14, 19, 26, 26, 29, 26, 9, 9, 14, 13, 19, 17, 15, 14, 20, 20, 3, 6, 1, 4, 31, 28, 1, 3, 2, 3, 3, 4, 11, 9, 4, 22, 23, 23, 21, 6, 7, 6, 26, 19, 24, 24, 28, 24, 14, 25, 22, 21, 21, 31, 16, 23, 21, 22, 21, 20, 32, 17, 21, 11, 13, 14, 24, 24, 22, 23, 24, 22, 25, 20, 19, 22, 26, 24, 20, 19, 28, 24, 23, 22, 12, 7, 11, 7, 16, 15, 4, 5, 1, 5, 11, 14, 27, 24, 24, 21, 21, 18, 22, 13, 7, 15, 4, 6, 12, 9, 10, 12, 8, 9, 25, 23, 24, 19, 20, 17, 30, 30, 27, 23, 28, 7, 3, 8, 9, 15, 31, 12, 14, 14, 17, 16, 18, 20, 18, 24, 21, 20, 22, 26, 21, 29, 5, 5, 5, 3, 8, 6, 9, 9, 5, 9, 14, 13, 8, 10, 10, 28, 28, 25, 27, 26, 9, 8, 8, 8, 6, 11, 30, 25, 28, 27, 27, 1, 3, 25, 23, 26, 25, 21, 22, 19, 19, 29, 17, 20, 13, 14, 13, 16, 8, 12, 9, 17, 14, 18, 13, 10, 22, 23, 23, 3, 2, 7, 1, 16, 17, 21, 17, 7, 6, 6, 2, 2, 8, 9, 17, 21, 17, 14, 19, 17, 23, 7, 23, 26, 30, 22, 18, 6, 4, 6, 4, 4, 6, 5, 30, 11, 10, 8, 9, 12, 12, 9, 8, 8, 11, 21, 14, 17, 27, 7, 14, 12, 14, 18, 11, 12, 9, 10, 9, 10, 7, 8, 18, 16, 18, 25, 25, 31, 7, 4, 7, 7, 10, 8, 14, 17, 21, 19, 29, 18, 18, 30, 29, 29, 18, 20, 10, 15, 10, 13, 14, 6, 7, 12, 9, 11, 26, 27, 31, 25, 25, 31, 28, 27, 30, 28, 13, 32, 28, 9, 12, 15, 13, 8, 17, 21, 15, 14, 15, 15, 27, 29, 27, 23, 14, 9, 16, 18, 28, 26, 29, 29, 26, 8, 26, 30, 12, 24, 26, 24, 4, 6, 23, 23, 17, 18, 22, 20, 25, 35, 11, 16, 11, 4, 4, 6, 6, 9, 4, 3, 9, 32, 26, 29, 24, 22, 8, 22, 12, 18, 16, 18, 18, 16, 2, 6, 22, 9, 20, 24, 20, 20, 21, 18, 4, 8, 2, 5, 5, 8, 15, 18, 15, 17, 13, 13, 8, 6, 12, 14, 16, 13, 27, 20, 24, 21, 26, 21, 20, 23, 18, 17, 14, 7, 8, 12, 8, 7, 16, 22, 22, 24, 19, 20, 25, 24, 21, 22, 22, 19, 20, 25, 27, 27, 24, 24, 23, 25, 23, 22, 23, 27, 25, 23, 27, 19, 19, 15, 19, 24, 24, 21, 25, 27, 24, 21, 26, 24, 25, 25, 24, 25, 23, 28, 35, 28, 28, 30, 8, 11, 7, 6, 9, 22, 19, 21, 24, 19, 24, 27, 25, 25, 13, 18, 12, 18, 5, 7, 8, 7, 22, 26, 27, 23, 24, 11, 12, 7, 27, 30, 31, 15, 7, 6, 4, 7, 27, 25, 25, 23, 25, 17, 7, 4, 18, 18, 16, 17, 18, 15, 18, 30, 30, 27, 28, 16, 11, 14, 28, 15, 12, 10, 14, 12, 28, 23, 28, 28, 7, 10, 26, 26, 26, 22, 25, 3, 12, 9, 14, 6, 9, 6, 22, 23, 11, 10, 13, 11, 13, 13, 11, 17, 11, 14, 13, 11, 12, 13, 22, 22, 22, 21, 35, 13, 21, 22, 20, 16, 18, 8, 7, 2, 6, 7, 9, 8, 9, 13, 14, 14, 17, 9, 11, 11, 5, 8, 8, 24, 21, 29, 22, 24, 28, 27, 23, 27, 23, 26, 28, 14, 13, 12, 7, 11, 21, 19, 28, 28, 28, 9, 7, 10, 7, 10, 9, 20, 22, 21, 19, 22, 22, 15, 20, 19, 23, 24, 21, 19, 22, 24, 24, 24, 28, 27, 25, 7, 7, 26, 25, 22, 35, 37, 35, 31, 11, 9, 10, 12, 19, 21, 19, 21, 7, 3, 5, 9, 11, 14, 13, 13, 15, 15, 13, 12, 21, 19, 20, 19, 19, 26, 26, 24, 22, 26, 3, 35, 19, 22, 21, 20, 20, 10, 6, 5, 6, 6, 5, 27, 25, 24, 1, 3, 14, 12, 10, 19, 16, 16, 16, 14, 4, 6, 6, 2, 7, 24, 22, 5, 9, 11, 8, 18, 20, 22, 24, 27, 29, 26, 27, 19, 16, 28, 32, 33, 18, 16, 23, 20, 22, 20, 18, 7, 6, 22, 32, 28, 30, 31, 6, 4, 4, 23, 23, 17, 17, 20, 19, 19, 32, 32, 1, 6, 3, 34, 33, 30, 6, 17, 17, 4, 5, 4, 24, 24, 33, 32, 31, 24, 25, 27, 24, 6, 8, 4, 27, 16, 18, 6, 9, 24, 22, 25, 23, 5, 4, 5, 12, 4, 4, 7, 7, 4, 5, 10, 11, 13, 12, 3, 8, 6, 7, 15, 13, 20, 5, 2, 19, 21, 24, 21, 38, 18, 21, 31, 32, 7, 4, 22, 19, 30, 28, 13, 25, 4, 17, 15, 15, 28, 8, 7, 8, 9]
    # fmt: on
    return (b,)


@app.cell
def _(b, np, plt):
    def _():
        fig, ax1 = plt.subplots(1, 1, figsize=(20, 5))
        ax: plt.Axes = ax1
        ax.set_xlabel("idx")
        ax.set_ylabel("Q-score")

        ax.plot(np.arange(len(b)), b, ls="--", color="grey", mfc="blue", mec="blue", marker=".")
        ax.set_ylim(0, 40)
        ax.set_xlim(0, 856)

        fig.tight_layout()
        return fig


    _()
    return


@app.cell
def _(pymsaviz):
    def _():
        msa_file = open("data/3.fasta", "r")
        mv = pymsaviz.MsaViz(msa_file, show_consensus=True, show_count=True)
        return mv.savefig("out2.png")


    _()
    return


@app.cell
def _(y8):
    y8
    return


if __name__ == "__main__":
    app.run()
