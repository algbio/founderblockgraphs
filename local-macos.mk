CLANG_ROOT	= /usr/local/homebrew/opt/llvm
SDSL_ROOT	= /Users/tsnorri/Checkout/sdsl-lite/build

CXX			= $(CLANG_ROOT)/bin/clang++
CPPFLAGS	= -I$(CLANG_ROOT)/include/c++/v1 -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk
