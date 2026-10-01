<div align="center">

# Introduction to Compressed Sensing

From signal transforms to sparse recovery, in one tutorial.

</div>

---

## About

Compressed sensing asks a question that sounds impossible at first: if a signal is sparse
in some basis, can you reconstruct it exactly from far fewer samples than the Nyquist rate
demands? The answer is yes, and this tutorial builds up to it from the ground.

The path is deliberate. It starts with signal transformations, because sparsity only means
something relative to a basis, and it covers the Fourier, wavelet and curvelet families. It
then introduces the concepts the theory needs: sparsity, norms, and convex optimization.
Only after that does it develop compressed sensing itself, along with the observation
matrices that make recovery possible and the greedy algorithms that carry it out. The final
chapter applies all of it to images.

The Chinese name for the subject, 压缩感知, was coined by Professor Qionghai Dai of
Tsinghua University.

## Who this is for

You should have taken a signals and systems course and a linear algebra course first. The
tutorial assumes you know the Fourier transform, and it uses matrix notation freely. No
prior exposure to compressed sensing is assumed.

## Contents

### 1. Signal Transformation

Sparsity only means something relative to a basis, so the tutorial starts by building the
vocabulary of transforms. Eight of them, in two families. The Fourier family covers the
discrete Fourier transform, its non-uniform variant, and the discrete cosine transform. The
time-frequency family covers the short time Fourier transform and the continuous and
discrete wavelet transforms. The chapter closes with the continuous and discrete curvelet
transforms, which handle edges better than wavelets do. Each is given its definition, a
fast algorithm, and a reason to choose it over the others.

### 2. Basic Concepts

The vocabulary the theory runs on. Sparsity and how it is measured; norms, including why
the l1 norm behaves so differently from the l0 "norm"; and convex optimization as the
machinery that turns a sparse-recovery idea into something solvable. The chapter closes
with three applied tools the later chapters depend on: projection matrices, image quality
assessment, and BayesShrink.

### 3. Compressed Sensing

The theory itself. The recovery condition, the restricted isometry property, and the
observation matrices that satisfy it in practice, including random Gaussian and Bernoulli
matrices. The last section covers greedy pursuit: orthogonal matching pursuit first, then
CoSaMP and SAMP, which correct earlier mistakes and need less prior knowledge.

### 4. Image processing

The application. Total variation regularization penalises how much an image changes from
one pixel to the next rather than how many pixels are non-zero, which preserves edges that
a sparsity penalty would smear. Then ADMM, the splitting method used to solve it.

### Afterwards and References

Closing thoughts, and the bibliography.

## What this tutorial does not cover

- **Convex optimization in depth.** Chapter 2 introduces what compressed sensing needs and
  no more. For the full treatment, see the companion *Introduction to Convex Optimization*.
- **Non-greedy recovery.** The reconstruction algorithms here are greedy pursuit methods.
  Linear-programming approaches such as basis pursuit are mentioned where they matter but
  are not developed.
- **Two-dimensional extension beyond images.** The final chapter applies the theory to
  images, not to general multi-dimensional or streaming settings.

## Building the PDF

You need a TeX distribution with XeLaTeX. A **full** TeX Live installation is required
rather than a minimal one. The document class is bundled in this repository, so there is
nothing extra to install. (Tested with TeX Live 2023.)

```bash
latexmk
```

A `.latexmkrc` is included, so `latexmk` selects XeLaTeX automatically and runs as many
passes as the cross-references need.

```bash
latexmk -c    # remove auxiliary files, keep the PDF
latexmk -C    # remove everything, including the PDF
```

The output is `Introduction to Compressed Sensing.pdf`.

## Repository layout

```
Introduction to Compressed Sensing.tex   the tutorial itself
elegantbook.cls                          document class
assets/                                  figures used in the text
.latexmkrc, .gitignore                   build configuration
```

## License

This tutorial is released under
[CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/). See [LICENSE](LICENSE).
