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
    return SeqIO, mo, np, pl, plt, pymsaviz, sns, stats


@app.cell
def _(pl):
    data_df = pl.read_csv("data.csv")
    data_df
    return (data_df,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Quality Score Histogram""")
    return


@app.cell(hide_code=True)
def _(SeqIO, data_df, mo, np, pl, plt):
    def _():
        df = data_df.filter(pl.col("filepath").str.ends_with(".fastq"))

        qs_counts = np.zeros(128)  # Phred quality score goes from 0 to 40
        file_qs_avg_seq_mean = []
        file_qs_lowest_seq_mean = []

        for idx, filepath in enumerate(df["filepath"]):
            filestem = df["filestem"][idx]
            file_qs_seq_means = []
            # fastq is equivalent to fastq-sanger → Phred quality score
            for record in SeqIO.parse(filepath, "fastq"):
                q_scores = np.array(record.letter_annotations["phred_quality"])
                file_qs_seq_means.append(np.mean(q_scores))
                qs_counts += np.bincount(q_scores, minlength=128)

            file_qs_avg_seq_mean.append(np.mean(file_qs_seq_means))
            file_qs_lowest_seq_mean.append(np.min(file_qs_seq_means))

        info_df = pl.DataFrame(
            data={
                "filestem": df["filestem"],
                "avg mean seq Q-score": file_qs_avg_seq_mean,
                "min mean seq Q-score": file_qs_lowest_seq_mean,
            }
        )

        freqs = qs_counts / np.sum(qs_counts)

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 5))
        fig.suptitle("Q-score distribution")

        ax: plt.Axes = ax1
        ax.set_xlabel("Q-score")
        ax.set_ylabel("frequency")
        ax.set_xlim(0, 70)

        ax.bar(np.arange(128), freqs)

        ax: plt.Axes = ax2
        ax.set_xlabel("Q-score (∈ {0,...,40} and normalized)")
        ax.set_ylabel("frequency")
        ax.set_xlim(0, 41)

        freqs_clean = freqs[0:41] / np.sum(freqs[0:41])
        ax.bar(np.arange(41), freqs_clean)

        return freqs_clean, info_df, fig


    Q_score_freqs, _info_df, _fig = _()
    Q_score_freqs_cummul = np.cumulative_sum(Q_score_freqs)
    print("## Cumulative sum of Q-score frequences (where Q-score ∈ {0,...,40} and normalized)")
    print(Q_score_freqs_cummul)

    mo.vstack([_fig, _info_df])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Read Length Histogram""")
    return


@app.cell(hide_code=True)
def _(SeqIO, data_df, np, pl, plt, sns):
    def _():
        df = data_df.filter(pl.col("filepath").str.ends_with(".fastq"))

        def gen_tdf(filepaths) -> pl.DataFrame:
            max_count = 1500
            rl_counts = np.zeros(max_count + 1)
            for filepath in filepaths:
                # fastq is equivalent to fastq-sanger → Phred quality score
                counts = np.bincount(
                    [
                        len(record.seq) if len(record.seq) < max_count else max_count
                        for record in SeqIO.parse(filepath, "fastq")
                    ],
                    minlength=max_count + 1,
                )
                rl_counts += counts

            freqs = rl_counts / np.sum(rl_counts)

            return pl.DataFrame(
                data={
                    "Read length (bp)": np.arange(max_count + 1),
                    "freq": freqs,
                }
            )

        fig, axes = plt.subplots(2, 2, figsize=(17, 7))

        ax: plt.Axes = axes[0][0]
        ax.set_title("All (Any combination)")
        tdf = gen_tdf(df["filepath"])
        sns.lineplot(data=tdf, x="Read length (bp)", y="freq", ax=ax)

        ax: plt.Axes = axes[0][1]
        ax.set_title("rbcL only")
        tdf = gen_tdf(
            df.filter(
                pl.col("matk") == False,
                pl.col("rbcL") == True,
                pl.col("psbA-trnH") == False,
                pl.col("ITS") == False,
            )["filepath"]
        )
        sns.lineplot(data=tdf, x="Read length (bp)", y="freq", ax=ax)

        ax: plt.Axes = axes[1][0]
        ax.set_title("psbA-trnH only")
        tdf = gen_tdf(
            df.filter(
                pl.col("matk") == False,
                pl.col("rbcL") == False,
                pl.col("psbA-trnH") == True,
                pl.col("ITS") == False,
            )["filepath"]
        )
        sns.lineplot(data=tdf, x="Read length (bp)", y="freq", ax=ax)

        ax: plt.Axes = axes[1][1]
        ax.set_title("ITS only")
        tdf = gen_tdf(
            df.filter(
                pl.col("matk") == False,
                pl.col("rbcL") == False,
                pl.col("psbA-trnH") == False,
                pl.col("ITS") == True,
            )["filepath"]
        )
        sns.lineplot(data=tdf, x="Read length (bp)", y="freq", ax=ax)

        fig.tight_layout()
        return fig


    _()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Read Length Fit""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    Let a discrete random variable $X$, distributed by a negative binomial, we can write $X \sim NBin(r, p)$

    This $(r, p)$ parametrization is quite common, if $r$ is an integer, then the negative binomial distribution can be thought of as the distribution of the number of failures in a sequence of Bernoulli trials that continue until `r` successes occur, where `p` is the probability of success in a single Bernoulli trial.

    With the above parametrization, the negative log-likelihood function is:  
    $$\ln L(r, p | y) = \sum_{i} \left\{ \ln \Gamma(y_i + r) - \ln \Gamma(r) - \ln \Gamma(y_i + 1) + r \ln(p) + y_i \ln(1 - p) \right\}$$

    We can also parametrize using the mean $\mu$, and the variance $\sigma^2 \:$:

    $$p=\frac{\mu}{\sigma^2} \quad r=\frac{\mu^2}{\sigma^2 - \mu}$$
    """
    )
    return


@app.cell(hide_code=True)
def _(SeqIO, data_df, np, pl, plt, stats):
    def _():
        from scipy.special import gammaln
        from numpy import log as ln
        from scipy.optimize import minimize

        def nbin_nll(params, counts):
            """
            Negative log-likelihood for negative binomial distribution
            Parameters
            ----------
            params : tuple
                tuple containing r and p
            counts : np.ndarray
                array of shape (N,) containing counts.

            Returns
            -------
            float
                Negative log-likelihood
            """
            r, p = params
            if r <= 0 or p <= 0 or p >= 1.0:  # Parameter constraints
                return np.inf
            x = counts
            return -np.sum(gammaln(x + r) - gammaln(r) - gammaln(x + 1) + r * ln(p) + x * ln(1 - p))

        df = data_df.filter(pl.col("filepath").str.ends_with(".fastq"))

        def fit_and_vis(filepaths, ax: plt.Axes):
            data = np.concat(
                [
                    np.array([len(record.seq) for record in SeqIO.parse(filepath, "fastq")])
                    for filepath in filepaths
                ],
                axis=0,
            )
            r_mle, p_mle = minimize(fun=nbin_nll, x0=[1, 0.5], args=(data,)).x

            x = np.arange(1500)
            data_bincount = np.bincount(data, minlength=1500)[:1500]
            y_act = data_bincount / np.sum(data_bincount)
            y_fit = stats.nbinom.pmf(x, r_mle, p_mle)

            ax.plot(x, y_act, label="Actual")
            ax.plot(x, y_fit, label=f"PMF of $NBin(r={r_mle:.3E}, p={p_mle:.3E})$")
            ax.legend(loc="upper right")

        fig, axes = plt.subplots(2, 2, figsize=(17, 7))

        ax: plt.Axes = axes[0][0]
        ax.set_title("All (Any combination)")
        fit_and_vis(df["filepath"], ax)

        ax: plt.Axes = axes[0][1]
        ax.set_title("rbcL only")
        fit_and_vis(
            df.filter(
                pl.col("matk") == False,
                pl.col("rbcL") == True,
                pl.col("psbA-trnH") == False,
                pl.col("ITS") == False,
            )["filepath"],
            ax,
        )

        ax: plt.Axes = axes[1][0]
        ax.set_title("psbA-trnH only")
        fit_and_vis(
            df.filter(
                pl.col("matk") == False,
                pl.col("rbcL") == False,
                pl.col("psbA-trnH") == True,
                pl.col("ITS") == False,
            )["filepath"],
            ax,
        )

        ax: plt.Axes = axes[1][1]
        ax.set_title("ITS only")
        fit_and_vis(
            df.filter(
                pl.col("matk") == False,
                pl.col("rbcL") == False,
                pl.col("psbA-trnH") == False,
                pl.col("ITS") == True,
            )["filepath"],
            ax,
        )

        fig.tight_layout()
        return fig


    _()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Q-score block visualization""")
    return


@app.cell(hide_code=True)
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

        fig.tight_layout()
        return fig


    _()
    return


@app.cell(hide_code=True)
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

        fig, ax1 = plt.subplots(1, 1, figsize=(8, 3))
        ax: plt.Axes = ax1
        ax.set_xlabel("frequency")
        ax.set_ylabel("block size")
        plt.gca().invert_yaxis()
        ax.barh(np.arange(30), (counts / np.sum(counts))[:30])

        freqs_safe = counts[:15] / np.sum(counts[:15])
        print(freqs_safe)
        freqs_safe_cummul = np.cumulative_sum(freqs_safe)
        print(freqs_safe_cummul)

        fig.tight_layout()
        return fig


    _()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Negative binomial distribution visualization""")
    return


@app.cell(hide_code=True)
def _(np, plt, stats):
    def vis_nbin(r, p):
        rv = stats.nbinom(r, p)
        fig, ax1 = plt.subplots(1, 1, figsize=(12, 4))

        ax: plt.Axes = ax1
        ax.set_title(f"$X \\sim NBin({r:.3E}, {p:.3E})$")
        x = np.arange(0, rv.ppf(0.99))
        y = rv.pmf(x)

        ax.vlines(x, 0, y, colors="k", linestyles="-", lw=1)
        ax.plot(x, y, "bo", ms=2)

        fig.tight_layout()
        return fig
    return (vis_nbin,)


@app.cell(hide_code=True)
def _(mo):
    r_slider = mo.ui.slider(start=0.5, stop=100, step=0.5, show_value=True, label="r")
    p_slider = mo.ui.slider(start=0.0, stop=1, step=0.001, show_value=True, label="p")

    r_slider, p_slider
    return p_slider, r_slider


@app.cell(hide_code=True)
def _(p_slider, r_slider, vis_nbin):
    def _():
        r, p = r_slider.value, p_slider.value
        return vis_nbin(r, p)


    _()
    return


@app.cell(hide_code=True)
def _(mo):
    mean_slider = mo.ui.slider(start=1, stop=1001, step=5, show_value=True, label="mean")
    std_slider = mo.ui.slider(start=1, stop=2001, step=5, show_value=True, label="std")

    mean_slider, std_slider
    return mean_slider, std_slider


@app.cell(hide_code=True)
def _(mean_slider, std_slider, vis_nbin):
    def _():
        mean, std = mean_slider.value, std_slider.value
        r = mean * mean / (std * std - mean)
        p = mean / (std * std)
        return vis_nbin(r, p)


    _()
    return


@app.cell
def _(mo):
    val_slider = mo.ui.slider(start=200, stop=2000, step=10, show_value=True, label="length")
    val_to_mean = 0.4
    val_to_std = 1

    val_slider
    return val_slider, val_to_mean, val_to_std


@app.cell(hide_code=True)
def _(val_slider, val_to_mean, val_to_std, vis_nbin):
    def _():
        val = val_slider.value
        mean = val_to_mean * val
        std = val_to_std * val
        r = mean * mean / (std * std - mean)
        p = mean / (std * std)
        return vis_nbin(r, p)

    _() 
    return


@app.cell(disabled=True)
def _(pymsaviz):
    def _():
        msa_file = open("data/3.fasta", "r")
        mv = pymsaviz.MsaViz(msa_file, show_consensus=True, show_count=True)
        return mv.savefig("out2.png")


    _()
    return


if __name__ == "__main__":
    app.run()
