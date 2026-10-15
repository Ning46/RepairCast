# RepairCast

**Predicting CRISPR-Cas9 repair outcomes from the guide sequence alone.**

When Cas9 cuts DNA, the cell repairs the double-strand break, mostly through non-homologous end joining, and leaves behind a spread of insertions and deletions (indels). That spread is not random: it depends heavily on the sequence around the cut site. RepairCast takes a 23-nt target sequence (20-nt protospacer + `NGG` PAM) and predicts four summary statistics of the resulting repair profile.

| Target | Meaning | Range in training data |
|---|---|---|
| `Fraction_Insertions` | Share of edited reads that carry an insertion | 0 – 1 |
| `Avg_Deletion_Length` | Mean length of deletions (bp) | ~2 – 46 |
| `Indel_Diversity` | Entropy-style diversity of distinct indel outcomes | ~1 – 5.5 |
| `Fraction_Frameshifts` | Share of outcomes that shift the reading frame | 0 – 1 |

These matter in practice: a guide whose repairs are mostly frameshifts is a better choice for a gene knockout, and a guide with low indel diversity gives more predictable edits.

## Approach

The full analysis is in [`crispr_repair_prediction.ipynb`](crispr_repair_prediction.ipynb).

**Data exploration.** 1,065 training and 304 test sequences, all exactly 23 nt, no missing values. The four targets are only weakly correlated (max |r| = 0.52), so each is modelled separately.

**Feature engineering.** Four sequence encodings were compared on a held-out validation split:

1. **3-mer counts**: 64 count features per sequence. Loses positional information.
2. **Position-wise one-hot encoding**: each of the 23 positions encoded as A/C/G/T. Positions 22–23 are always `G` (the PAM), so 86 informative columns remain.
3. **Word2Vec 3-mer embeddings** (gensim, 50-d, averaged). Little to learn from such short sequences, and it performed worst.
4. **One-hot + GC content**: approach 2 plus the sequence's GC fraction. **Best, and used for the final model.**

One-hot encoding works well here because repair outcomes depend on *which base sits at which position* relative to the cut site (3 bp upstream of the PAM). Microhomology around the break, for example, strongly drives deletion length. Count- and embedding-based features discard exactly that positional signal.

**Model.** `MultiOutputRegressor(GradientBoostingRegressor)`, one gradient-boosted tree ensemble per target on standardised features. Random forest, AdaBoost, SVR, Bayesian ridge and XGBoost were also tried; gradient boosting scored best on the competition leaderboard.

## Results

80/20 train-validation split (`random_state=42`):

| Model | Validation R² | Validation MSE |
|---|---|---|
| GradientBoosting, default parameters | **0.227** | 2.50 |
| Best grid point (`n_estimators=100`, `lr=0.05`) | 0.227 | – |

A grid search over `n_estimators ∈ {100, 200, 300, 500}` × `learning_rate ∈ {0.01, 0.05, 0.1}` gave validation R² between 0.13 and 0.23. The tuned settings did not beat the defaults on the online leaderboard, so the final model keeps the default parameters and is retrained on the full training set.

> MSE is averaged over the four targets on their raw scales, so it is dominated by `Avg_Deletion_Length`. R² is the more comparable number.

## Run it

```bash
git clone https://github.com/Ning46/RepairCast.git
cd RepairCast
pip install -r requirements.txt
```

The competition data is **not included** in this repository. Competition page is not there, which is a private competition organized by company:

```bash
jupyter notebook crispr_repair_prediction.ipynb
```

Running all cells writes `submission.csv` with predictions for the test set.

## Limitations and next steps

- **Small data.** About 1,000 guides is little for a 4-target regression; R² ≈ 0.23 shows that most of the variance is still unexplained.
- **Single split.** Model selection used one 80/20 split, which is noisy at this size. K-fold cross-validation would give steadier estimates and explain why the validation and leaderboard rankings disagreed.
- **Biologically informed features.** Explicit microhomology features (length and position of repeated sequence flanking the cut), melting temperature and dinucleotide one-hot encodings are known to predict indel profiles and are the natural next step. Published models such as inDelphi and FORECasT are built around them.
- **Robust test encoding.** Test one-hot columns are currently built independently of training. They align on this dataset, but reindexing to the training columns would make the pipeline safe on new sequences.

## License

MIT
