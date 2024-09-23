# libertiff

## What is libertiff?

Libertiff is a C++11 simple, [header-only](libertiff.hpp), TIFF reader. It is MIT licensed.

Handles both ClassicTIFF and BigTIFF, little-endian or big-endian ordered.

The library does not (yet?) offer codec facilities. It is mostly aimed at
browsing through the linked chain of Image File Descriptors and their tags.

The library is thread-safe (that is the instances that it returns can
be used from multiple threads), if passed FileReader instances are themselves
thread-safe.

The library does not throw exceptions (but underlying std library might
throw exceptions in case of out-of-memory when allocating the tag array)

Optional features:
- define LIBERTIFF_C_FILE_READER before including that file, so that
  the libertiff::CFileReader class is available

## How to use it?

Look at the [demo.cpp](demo.cpp) test program.
