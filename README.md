# atomic-halide package

An Atom editor package for rapid development of simple Halide imaging language pipelines. See [this demo][demo] to see it in action.

Purpose
-------
This package was the direct result of wanting a more streamlined way to learn the basics of image processing using [Halide][halide]. Halide has an excellent set of well commented tutorials, but as I was experimenting with them, I found myself wanting a tighter iteration process than running make and looking at written image files after each code change.

This is noted as much as anything to say that this is a prototype of a live coding environment for imaging algorithms and not a finished product. It was created for the purposes of learning a number of technologies (Atom's Electron app framework, Node.js, image processing in general, Halide in particular, using web technology for desktop app UI, CoffeeScript, etc). I knew very little about any of these when I started. The UI is also very much of the "simplest thing that could possibly work" variety. I am neither a UI or web developer.

All that is to say that if there are rookie mistakes in the code or examples, feel free to point them out (politely).

Installation/Setup
------------------

Prerequisites:
* A copy of [Github's Atom Editor][atom] (most recently tested with 1.23.3).
* A C++ build toolchain: Xcode on the Mac, Visual Studio on Windows, make and g++ on Linux.
	With VS2017 you need to install the v140 toolset as well and use a "Developer command prompt for VS2017".
* A [binary distribution of Halide][halide-release] or a compatible version [built from source][halide-git].

	Note that because Atomic Halide dynamically links to compiled Halide filters and calls them with an FFI, you must pick a binary distribution whose architecture that matches that of Atom itself. As of this writing, that is 64-bit on Mac, Linux and Windows.
* A working version of Python (used by node-gyp to install native node modules via Atom's package manager).

	I used 32-bit Python 2.7.10 from [python.org][python] on Windows and the pre-installed versions of Python on Mac and Linux). As of this writing, node-gyp is not compatible with Python 3 and its error message when failing to build due to a missing or incorrect version of Python is not obvious.
* Clone this git repo.

    On Mac and Linux, you can put it anywhere and run setup.sh to symlink it in as a local Atom package and run apm install to fetch and build the native modules.

    On Windows, I found it is easier to just put the repo in Atom's local package directory directly (.atom/packages in your home directory) and run setup.bat from there to perform the apm install.
* Restart Atom or run the Reload command in the View menu to pick up the newly installed package.
* Add the examples folder of the git repo as a project folder.
* Open one of the cpp files (stretch.cpp, mosiac.cpp, etc).
* In the Packages menu for Atomic Halide, invoke the Toggle Preview command.
* On first use, press the Configure button at the top of the panel and select the folder where your build of Halide resides (it should have the include and bin or Debug/Release subfolders).

How it Works
------------
When the preview is opened, a new editor is opened, or the current editor is saved, the package attempts to shell out to make (Mac/Linux) or nmake (Windows) to build a shared library (dylib, so, dll on Mac, Linux, and Windows respectively). The details of how to go from a source file to a dylib are mostly defined by target definitions in the Makefile/NMakefile in the examples folder.

If there is no target, it leaves the existing panel in place. If there is a target and it fails, the error output is presented in the panel. If it succeeds, a JavaScript FFI is used to bind to the render function and parse the exported parameter metadata to build the panel of controls (sliders, checkboxes, and text fields).

As those inputs are manipulated, the render function is invoked (in process) and the buffer is copied into an ImageData to paint the result onto the HTML5 canvas element. If the source file is changed and saves, the panel is regenerated and the inputs are reset.

This really streamlined my workflow for playing with imaging algorithms. It makes it easy to introduce new parameters and see how various values affect the output. That iterative feedback loop makes the process much more fun.

Thanks
------
Special thanks go out to those who've contributed to Halide itself. It's a really neat piece of work and the excellent tutorials made it quite approachable. Thanks also to Github for Atom and the various authors of the [Node modules I used](package.json). Having such excellent building blocks made this much easier to assemble.

Caveats
-------
Most of the work on this package was done on a Mac. I did get it working on Windows as well, but I haven't exposed the required toolchain paths as configuration, so a different version of Visual Studio and the Windows SDK probably requires a slight tweak to a few constants in the [package source code](lib/winenv.coffee) right now. It has also been tested on Ubuntu 14.04 and works there.

Also, the package does nothing to set up the build toolchain (make, a C++ compiler, etc). It assumes it can shell out to make or nmake and that those in turn can determine the correct compiler to use.

Since the FFI runs the render function synchronously in process, you can hang or crash Atom if you make serious mistakes in your imaging code. There are certain kinds of bound checking mistakes that Halide will report via an error callback and those are handled more gracefully. I experimented with various degrees of async rendering (out of process, in process but using the FFI async call path), but haven't made the switch yet because I didn't like the way they felt to use, so I haven't done that yet.

These and a few others are also noted in the [TODO list](TODO.txt).

[atom]: https://atom.io
[demo]: https://www.youtube.com/watch?v=PTSVlT3Iq4U
[halide]: http://halide-lang.org/
[halide-git]: https://github.com/halide/Halide
[halide-release]: https://github.com/halide/Halide/releases
[python]: https://www.python.org
