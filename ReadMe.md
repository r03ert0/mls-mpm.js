# A Javascript version of Taichi's 88-line material point method

Roberto Toro, March 2019

<img src="https://user-images.githubusercontent.com/2310732/53946632-7d367780-40c4-11e9-8ceb-93772240351b.gif" width="400px">

This is a javascript version of the awesome implementation in 88 lines of c++ of the Moving Least Squares Material Point Method, by Yuanming Hu.

The original repository is https://github.com/yuanming-hu/taichi_mpm, and the c++ code can be directly downloaded from http://bit.ly/mls-mpm88. This link contains 2 files: `taichi.h`, with linear algebra and helper functions; and `mls-mpm88.cpp` with the actual implementation of the material point method.

In this javascript version, the linear algebra functions are in the `algebra.js` file, and the MPM code in the `mls-mpm.js` file (with the `88` removed, because now it's 148 lines... but there's much more comments, so we're still good).

A main tricky thing in the translation, which still is not completely solved, is that matrices in Taichi are coded transposed. In the translation, there's still a mix of functions that work transposed (svd, polar_decomp, mulMatVec), and those that don't. That makes some parts of the code a bit convoluted. Ideally, I'd like to move all the code to standard non-transposed matrix encoding.

### Acknowledgements

* Thanks to [Yuanming Hu](https://github.com/yuanming-hu) for making the C++ code available and for help on understanding it.
* Thanks to [Kevin Chapelier](https://github.com/kchapelier) for making the javascript code x5 faster and reducing the initialisation from several seconds to barely noticeable.
