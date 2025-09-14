# myrio-core

### Computing results

let $s = (s_1, ..., s_n)$ the similarity scores between the cluster centroid and the k-mer counts of the entries

Softmax pooling applied to leaves:

$$
score = \frac{\sum_{i=1}^{n} s_i e^{\lambda s_i}}{\sum_{i=1}^{n} e^{\lambda s_i}} \cdot \frac{n}{n + \mu} \quad \text{where } \lambda = \lambda_\text{leaf}
$$

Softmax pooling applied to branches:

$$
\text{score} = \frac{\sum_{i=1}^{n} s_i e^{\lambda s_i}}{\sum_{i=1}^{n} e^{\lambda s_i}} \quad \text{where } \lambda = \lambda_\text{branch}
$$

Confidence score:

let $x_\text{first}$, $x_\text{second}$ the best and second best softmax pool scores at a given rank, let $\sigma_\text{pool}$ the standard deviation of the pooling scores

$$
x = \frac{x_\text{first} - x_\text{second}}{\sigma_\text{pool}}
$$
$$
\text{conf} =  x_\text{first}^\epsilon \cdot \frac{x^\gamma}{x^\gamma + \delta}
$$
