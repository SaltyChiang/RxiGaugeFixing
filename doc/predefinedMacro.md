# Macros to make sure how to perform

## Checking for OS (platform)

    Linux and Linux-derived           __linux__
    Android                           __ANDROID__ (implies __linux__)
    Linux (non-Android)               __linux__ && !__ANDROID__
    Darwin (Mac OS X and iOS)         __APPLE__
    Akaros (http://akaros.org)        __ros__
    Windows                           _WIN32
    Windows 64 bit                    _WIN64 (implies _WIN32)
    NaCL                              __native_client__
    AsmJS                             __asmjs__
    Fuschia                           __Fuchsia__

## Checking the compiler

    Visual Studio       _MSC_VER
    gcc                 __GNUC__
    clang               __clang__
    emscripten          __EMSCRIPTEN__ (for asm.js and webassembly)
    MinGW 32            __MINGW32__
    MinGW-w64 32bit     __MINGW32__
    MinGW-w64 64bit     __MINGW64__

## Checking compiler version

### gcc

`__GNUC__` (e.g. 5) and `__GNUC_MINOR__` (e.g. 1).

### clang

`__clang_major__`, `__clang_minor__`, `__clang_patchlevel__`.

### Visual Studio

`_MSC_VER` and `_MSC_FULL_VER`.

    VS                          _MSC_VER    _MSC_FULL_VER
    1                           800
    3                           900
    4                           1000
    4                           1020
    5                           1100
    6                           1200
    6   SP6                     1200        12008804
    7                           1300        13009466
    7.1 (2003)                  1310        13103077
    8   (2005)                  1400        140050727
    9   (2008)                  1500        150021022
    9   SP1                     1500        150030729
    10  (2010)                  1600        160030319
    10  (2010) SP1              1600        160040219
    11  (2012)                  1700        170050727
    12  (2013)                  1800        180021005
    14  (2015)                  1900        190023026
    14  (2015 Update 1)         1900        190023506
    14  (2015 Update 2)         1900        190023918
    14  (2015 Update 3)         1900        190024210
    15  (2017 Update 1 & 2)     1910        191025017
    15  (2017 Update 3 & 4)     1911
    15  (2017 Update 5)         1912

### MinGW

MinGW (aka MinGW32) and MinGW-w64 32bit: `__MINGW32_MAJOR_VERSION` and `__MINGW32_MINOR_VERSION`.

MinGW-w64 64bit: `__MINGW64_VERSION_MAJOR` and `__MINGW64_VERSION_MINOR`.

## Checking processor architecture

### gcc

- `__i386__`
- `__x86_64__`
- `__arm__`. If defined, you can further check:
  - `__ARM_ARCH_5T__`
  - `__ARM_ARCH_7A__`
- `__powerpc64__`
- `__aarch64__`