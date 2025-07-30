<div align="center">

# Introduction to Compressed Sensing
<hr>

</div>


This repository stores Compressed Sensing tutorial.

This tutorial uses [CC-BY-SA-4.0](https://creativecommons.org/licenses/by-sa/4.0/) License.

Human beings have run into 21st century. The rapid development of information technology
has greatly increased the amount of data we need. It's clear that the real world is analog
and the digital world is discrete so that signal sampling is the necessary way to convert
the real world into the digital one. **Shannon-Nyquist sampling theorem** tells us
that the sampling rate should be at least twice the maximum frequency of the signal to
avoid aliasing, and this frequency is called the Nyquist frequency. However,
Shannon-Nyquist sampling theorem is faced with a problem that us humans wants to decrease
the cost of data sampling and storage. Is there a way to break through the
Shannon-Nyquist sampling theorem and sample the signal at a lower rate? The answer is yes,
and this is the **Compressed Sensing**.

We often call Compressed Sensing as CS for short, and also Compressive Sampling. Its
Chinese name is **压缩感知**, translated by Professor Qionghai Dai from Tsinghua
University. Compressed Sensing is a new theory that can reconstruct the original signal
from a small number of samples, which is much less than the Nyquist rate. For example, if
we lost **70%** of the samples, we can still reconstruct the original signal with a
high probability. This is much unbelievable for us. Therefore someone said that Compressed
Sensing is the most important discovery in information theory since Shannon-Nyquist
sampling theorem.

To learn such a powerful theory is not easy. We need to start from the signal
transformation, and then some basic concepts like sparse representation and norm. The
final part is Compressed Sensing theory, we will introduce its concept and some
more like the observation matrices and the reconstruction algorithms. In this tutorial
we will mainly introduce greedy algorithms.

If you find any mistakes or have any suggestions, please feel free to contact me or
create an issue or a PR.
