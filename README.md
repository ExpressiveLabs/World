# TorchWorld - a Torch/C++ port and wrapper for [mmorise/World](https://github.com/mmorise/World)

TorchWorld is a Torch/C++ port and wrapper for the popular vocoder World by Masanori Morise.<br>**This repository is still a work-in-progress.** Don't expect wonders from it for now.

Besides rewriting functions to be compatible with Torch, this port also tries to follow modern C++ conventions where possible, for example using `std::shared_ptr` instead of raw pointers and `std::vector` instead of C-style arrays.

All rewritten functions and classes follow the same structure: they are contained within the `tw` namespace and have the same name as their original counterparts. Each class has a `run()` method that will execute its functionality. For example, to run CheapTrick, you would do the following:

```cpp
using namespace tw;

std::vector<double> x;
int fs = 44100;

auto options = std::make_shared<CheapTrick::Options>(fs);

CheapTrick cheaptrick(x, x.size(), fs, options);
cheaptrick.run();
```